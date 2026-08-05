/**
 * Locus Splice Graph Builder
 *
 * Transforms an Ensembl Gene (with Transcript/Exon via expand=1) into a
 * compact splice graph and supplies helpers to assemble transcript sequences.
 */
(function(global){
  function sortExonsForTranscript(tr) {
    const exons = (tr.Exon || []).slice();
    exons.sort((a,b) => tr.strand === -1 ? b.start - a.start : a.start - b.start);
    return exons;
  }

  function uniqueExonKey(e, strand) { return `${e.id || `${e.start}-${e.end}`}:${strand}`; }

  function buildSpliceGraph(gene) {
    const nodes = new Map(); // key -> {id,key,start,end,strand, transcripts:Set}
    const edges = new Map(); // `${fromKey}->${toKey}` -> {from,to, transcripts:Set}

    const transcripts = (gene.Transcript || []).map(tr => ({
      id: tr.id,
      name: tr.display_name || tr.id,
      strand: tr.strand || gene.strand,
      biotype: tr.biotype,
      is_canonical: !!tr.is_canonical,
      exons: sortExonsForTranscript(tr)
    }));

    for (const tr of transcripts) {
      const ordered = tr.exons;
      for (let i=0; i<ordered.length; i++) {
        const ex = ordered[i];
        const key = uniqueExonKey(ex, tr.strand);
        if (!nodes.has(key)) nodes.set(key, { id: ex.id || key, key, start: ex.start, end: ex.end, strand: tr.strand, transcripts: new Set() });
        nodes.get(key).transcripts.add(tr.id);
        if (i < ordered.length - 1) {
          const ex2 = ordered[i+1];
          const k2 = uniqueExonKey(ex2, tr.strand);
          const edgeKey = `${key}->${k2}`;
          if (!edges.has(edgeKey)) edges.set(edgeKey, { from: key, to: k2, transcripts: new Set() });
          edges.get(edgeKey).transcripts.add(tr.id);
        }
      }
    }

    // Serialize
    const nodeList = Array.from(nodes.values()).map(n => ({
      id: n.id,
      start: n.start,
      end: n.end,
      strand: n.strand,
      transcripts: Array.from(n.transcripts)
    }));
    const edgeList = Array.from(edges.values()).map(e => ({
      from: e.from,
      to: e.to,
      transcripts: Array.from(e.transcripts)
    }));

    // Identify start and end exons from annotated transcripts
    const startExonKeys = new Set();
    const endExonKeys = new Set();
    for (const tr of transcripts) {
      if (tr.exons.length > 0) {
        startExonKeys.add(uniqueExonKey(tr.exons[0], tr.strand));
        startExonKeys.add(tr.exons[0].id || uniqueExonKey(tr.exons[0], tr.strand));
        const last = tr.exons[tr.exons.length - 1];
        endExonKeys.add(uniqueExonKey(last, tr.strand));
        endExonKeys.add(last.id || uniqueExonKey(last, tr.strand));
      }
    }

    return { nodes: nodeList, edges: edgeList, transcripts, startExonKeys: Array.from(startExonKeys), endExonKeys: Array.from(endExonKeys) };
  }

  // Assemble DNA sequence for a transcript using a locus {dna,start,end,strand}
  function assembleTranscriptDNA(locus, transcript) {
    const parts = [];
    const exons = (transcript.exons || []).slice();
    for (const ex of exons) {
      // locus is oriented to genome +; exonSequenceFromLocus will reverse-complement if needed
      const seq = (global.LocusEnsembl && typeof global.LocusEnsembl.exonSequenceFromLocus === 'function')
        ? global.LocusEnsembl.exonSequenceFromLocus(locus, ex)
        : '';
      parts.push(seq);
    }
    return parts.join('');
  }

  function dnaToRNA(s) { return (s||'').toUpperCase().replace(/[^ACGT]/g,'').replace(/T/g,'U'); }

  function revComp(dna) {
    const comp = { A:'T', T:'A', C:'G', G:'C', N:'N' };
    return (dna||'').split('').reverse().map(b => comp[b] || 'N').join('');
  }

  function genomicSlice(locus, gStart, gEnd) {
    const start = Math.max(locus.start, gStart);
    const end = Math.min(locus.end, gEnd);
    if (end < start) return '';
    const sOff = start - locus.start;
    const eOff = end - locus.start;
    return locus.dna.substring(sOff, eOff + 1);
  }

  function orientedSeq(locus, gStart, gEnd) {
    const seq = genomicSlice(locus, gStart, gEnd);
    if (locus.strand === -1) return revComp(seq);
    return seq;
  }

  // Check canonical GT-AG signals at donor (after exon) and acceptor (before exon) in transcript orientation
  function hasCanonicalDonor(locus, exon) {
    if (locus.strand === 1) {
      const s = orientedSeq(locus, exon.end + 1, exon.end + 2);
      return s.length === 2 && s === 'GT';
    } else {
      const s = orientedSeq(locus, exon.start - 2, exon.start - 1);
      return s.length === 2 && s === 'GT';
    }
  }

  function hasCanonicalAcceptor(locus, exon) {
    if (locus.strand === 1) {
      const s = orientedSeq(locus, exon.start - 2, exon.start - 1);
      return s.length === 2 && s === 'AG';
    } else {
      const s = orientedSeq(locus, exon.end + 1, exon.end + 2);
      return s.length === 2 && s === 'AG';
    }
  }

  function deriveCombinatorialEdges(model, locus, rules = {}) {
    const { nodes } = model;
    const allowSkipEdges = rules.allowSkipEdges !== false; // default true
    const requireCanonical = rules.requireCanonicalSignals !== false; // default true

    const byPos = nodes.slice().sort((a,b)=> locus.strand===-1 ? b.start - a.start : a.start - b.start);
    const nodeKey = (n) => `${n.id || (n.start+'-'+n.end)}:${n.strand}`;
    const existing = new Set(model.edges.map(e=> `${e.from}->${e.to}`));
    const extra = [];

    // Annotated edges already in model.edges; add skip edges if allowed
    if (allowSkipEdges) {
      for (let i=0; i<byPos.length; i++) {
        for (let j=i+1; j<byPos.length; j++) {
          const a = byPos[i], b = byPos[j];
          // Ensure same strand and monotonic order in transcript direction
          if (a.strand !== b.strand) continue;
          const key = `${nodeKey(a)}->${nodeKey(b)}`;
          if (existing.has(key)) continue; // already annotated
          // Skip if directly overlapping in genomic coords (shouldn't happen for distinct exons)
          if (!(a.end < b.start || b.end < a.start)) continue;

          if (requireCanonical) {
            if (!hasCanonicalDonor(locus, a)) continue;
            if (!hasCanonicalAcceptor(locus, b)) continue;
          }

          extra.push({ from: nodeKey(a), to: nodeKey(b), transcripts: [] });
          existing.add(key);
        }
      }
    }
    return extra;
  }

  function buildEventModel(gene, locus, rules = {}) {
    const g = buildSpliceGraph(gene);
    const extra = deriveCombinatorialEdges(g, locus, rules);
    const edges = g.edges.concat(extra);
    return { ...g, edges };
  }

  function enumerateChains(eventModel, rules = {}) {
    const maxChains = rules.maxChains || 200;
    const maxExons = rules.maxChainExons || 12;
    const startSet = new Set(eventModel.startExonKeys);
    const endSet = new Set(eventModel.endExonKeys);

    // Build adjacency
    const adj = new Map();
    for (const e of eventModel.edges) {
      if (!adj.has(e.from)) adj.set(e.from, []);
      adj.get(e.from).push(e.to);
    }

    const chains = [];
    const visiting = new Set();

    function dfs(nodeKeyStr, path) {
      if (chains.length >= maxChains) return;
      if (path.length > maxExons) return;
      if (endSet.has(nodeKeyStr)) {
        chains.push(path.slice());
        // continue exploring to allow longer paths that also end at end exons? stop here to reduce explosion
      }
      const nbrs = adj.get(nodeKeyStr) || [];
      for (const nxt of nbrs) {
        const step = `${nodeKeyStr}->${nxt}`;
        if (path.includes(nxt)) continue; // avoid cycles
        path.push(nxt);
        dfs(nxt, path);
        path.pop();
        if (chains.length >= maxChains) return;
      }
    }

    // Seed DFS from starts
    const starts = Array.from(startSet);
    for (const s of starts) {
      dfs(s, [s]);
      if (chains.length >= maxChains) break;
    }

    // Convert keys back to exon objects with coords for assembly
    const keyToNode = new Map();
    for (const n of eventModel.nodes) {
      const k = `${n.id || (n.start+'-'+n.end)}:${n.strand}`;
      keyToNode.set(k, n);
    }
    const assembled = chains.map((route, idx) => ({ id: `CHAIN_${idx+1}`, route, exons: route.map(k=> keyToNode.get(k)).filter(Boolean) }));
    return assembled;
  }

  global.LocusSpliceGraph = {
    buildSpliceGraph,
    buildEventModel,
    enumerateChains,
    assembleTranscriptDNA,
    dnaToRNA
  };
})(typeof window !== 'undefined' ? window : globalThis);
