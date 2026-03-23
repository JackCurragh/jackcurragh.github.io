/**
 * SpliceModel - segments, edges, transcripts, collapsed mapping
 *
 * Converts Ensembl gene + locus into minimal exonic segments that reflect
 * truly shared nucleotides across transcripts. Builds adjacency among
 * segments (annotated + plausible edges) and per-transcript segment paths.
 */
(function(global){
  function revComp(dna){ const c={A:'T',T:'A',C:'G',G:'C',N:'N'}; return (dna||'').split('').reverse().map(b=>c[b]||'N').join(''); }
  function orientedSeq(locus, gStart, gEnd){
    const s = Math.max(locus.start, gStart), e = Math.min(locus.end, gEnd);
    if (e < s) return '';
    const so = s - locus.start, eo = e - locus.start;
    const slice = (locus.dna || '').substring(so, eo+1);
    return locus.strand === -1 ? revComp(slice) : slice;
  }
  function donorIsGT(locus, exon){
    if (locus.strand === 1){ return orientedSeq(locus, exon.end+1, exon.end+2) === 'GT'; }
    return orientedSeq(locus, exon.start-2, exon.start-1) === 'GT';
  }
  function acceptorIsAG(locus, exon){
    if (locus.strand === 1){ return orientedSeq(locus, exon.start-2, exon.start-1) === 'AG'; }
    return orientedSeq(locus, exon.end+1, exon.end+2) === 'AG';
  }

  function buildBreakpoints(gene){
    const bp = new Set();
    (gene.Transcript||[]).forEach(tr => {
      (tr.Exon||[]).forEach(ex => { bp.add(ex.start); bp.add(ex.end); });
    });
    return Array.from(bp).sort((a,b)=>a-b);
  }

  function buildSegments(gene, breakpoints){
    // Produce atomic segments by splitting every exon on global breakpoints.
    // Deduplicate by (start,end,strand).
    const key = (s,e,strand)=>`${s}-${e}:${strand}`;
    const segMap = new Map();
    (gene.Transcript||[]).forEach(tr => {
      const strand = tr.strand || gene.strand;
      const bps = breakpoints;
      (tr.Exon||[]).forEach(ex => {
        const s = Math.min(ex.start, ex.end), e = Math.max(ex.start, ex.end);
        for (let i=0;i<bps.length-1;i++){
          const a=bps[i], b=bps[i+1];
          // segment [a,b] is within exon if [a,b] subset [s,e]
          if (a < s || b > e) continue;
          const k = key(a,b,strand);
          if (!segMap.has(k)) segMap.set(k, { id: `SEG_${segMap.size+1}`, start: a, end: b, strand, sharedBy: new Set() });
          segMap.get(k).sharedBy.add(tr.id);
        }
      });
    });
    // Serialize
    const segments = Array.from(segMap.values()).map(s => ({ id: s.id, start: s.start, end: s.end, strand: s.strand, sharedBy: Array.from(s.sharedBy), label: `${s.start}-${s.end}` }));
    const byKey = new Map();
    segments.forEach(s => byKey.set(`${s.start}-${s.end}:${s.strand}`, s));
    const index = new Map(segments.map(s => [s.id, s]));
    return { segments, segmentIndex: index, byKey };
  }

  function transcriptPaths(gene, byKey){
    const paths = [];
    (gene.Transcript||[]).forEach(tr => {
      const strand = tr.strand || gene.strand;
      const exons = (tr.Exon||[]).slice().sort((a,b)=> strand===-1 ? b.start-a.start : a.start-b.start);
      const segIds = [];
      exons.forEach(ex => {
        const s = Math.min(ex.start, ex.end), e = Math.max(ex.start, ex.end);
        // walk across [s,e) by searching in byKey. As we split at all breakpoints, segments exist exactly.
        let cur = s;
        while (cur < e){
          const k = `${cur}-${Math.min(e, cur + (e - cur))}:${strand}`; // exact step unknown; find any segment starting at cur
          // Find segment starting at cur by brute scan of keys in byKey
          let seg = null;
          // Attempt exact hit first
          seg = byKey.get(`${cur}-${e}:${strand}`);
          if (!seg){
            // Try find next breakpoint end
            for (const [k2, s2] of byKey){
              const [range, st] = k2.split(':');
              if (parseInt(st,10) !== strand) continue;
              const [a,b] = range.split('-').map(Number);
              if (a === cur && b <= e){ seg = s2; break; }
            }
          }
          if (!seg) break; // defensive
          segIds.push(seg.id);
          cur = seg.end;
        }
      });
      if (segIds.length > 0){ paths.push({ id: tr.id, name: tr.display_name || tr.id, strand, segments: segIds }); }
    });
    return paths;
  }

  function buildAnnotatedEdges(paths){
    const edges = new Map();
    function addEdge(a,b,annotated){
      const k = `${a}->${b}`; if (!edges.has(k)) edges.set(k, { from:a, to:b, annotated: !!annotated, support: 0 });
      if (annotated) edges.get(k).annotated = true;
    }
    paths.forEach(p => {
      for (let i=0;i<p.segments.length-1;i++) addEdge(p.segments[i], p.segments[i+1], true);
    });
    return edges;
  }

  function addPlausibleEdges(edges, model, locus, rules){
    const { segments, segmentIndex } = model;
    if (!rules.allowSkipEdges) return;
    // Add potential skip edges from any segment to a downstream segment if canonical signals present at the exon boundaries.
    // Approximation: treat each segment as exon fragment; donor at segment end; acceptor at target segment start.
    // Require canonical GT-AG if rules.requireCanonicalSignals.
    const list = segments.slice().sort((a,b)=> a.start-b.start);
    for (let i=0;i<list.length;i++){
      const a = list[i];
      for (let j=i+1;j<list.length;j++){
        const b = list[j];
        if (a.strand !== b.strand) continue;
        if (a.end >= b.start) continue; // must be downstream and non-overlapping
        if (rules.requireCanonicalSignals){
          const donorOk = donorIsGT(locus, { start:a.start, end:a.end });
          const acceptorOk = acceptorIsAG(locus, { start:b.start, end:b.end });
          if (!donorOk || !acceptorOk) continue;
        }
        const k = `${a.id}->${b.id}`;
        if (!edges.has(k)) edges.set(k, { from:a.id, to:b.id, annotated:false, support:0 });
      }
    }
  }

  function unionExonicPieces(segments){
    const ints = segments.slice().sort((a,b)=> a.start-b.start).map(s => [s.start, s.end]);
    const pieces = [];
    for (const [s,e] of ints){
      if (!pieces.length || s > pieces[pieces.length-1][1]) pieces.push([s,e]);
      else pieces[pieces.length-1][1] = Math.max(pieces[pieces.length-1][1], e);
    }
    // annotate collapsed coords
    let c=0; const out=[];
    for (const [s,e] of pieces){ const len = Math.max(0, e-s); out.push({ start:s, end:e, cStart:c, cEnd:c+len }); c += len; }
    return out;
  }

  function g2cFactory(pieces){
    return function(g){
      // find piece containing g
      for (const p of pieces){ if (g >= p.start && g <= p.end) return p.cStart + (g - p.start); }
      // in intron: snap to nearest boundary
      let best=null, bestDist=Infinity, bestC=null;
      for (const p of pieces){
        if (g < p.start){ const d = p.start - g; if (d < bestDist){ bestDist=d; bestC=p.cStart; } }
        else if (g > p.end){ const d = g - p.end; if (d < bestDist){ bestDist=d; bestC=p.cEnd; } }
      }
      return bestC==null?0:bestC;
    };
  }

  function c2gFactory(pieces){
    return function(c){
      for (const p of pieces){ if (c >= p.cStart && c <= p.cEnd) return p.start + (c - p.cStart); }
      // outside: clamp
      if (c <= 0) return pieces[0]?.start || 0;
      return pieces[pieces.length-1]?.end || 0;
    };
  }

  function assignLanes(segments){
    const segs = segments.slice().sort((a,b)=> a.start - b.start);
    const lanes = []; // each lane has lastEnd
    const laneOf = new Map();
    for (const s of segs){
      let placed = false;
      for (let i=0;i<lanes.length;i++){
        if (s.start > lanes[i]){ laneOf.set(s.id, i); lanes[i] = s.end; placed=true; break; }
      }
      if (!placed){ laneOf.set(s.id, lanes.length); lanes.push(s.end); }
    }
    return { laneOf, laneCount: lanes.length };
  }

  function buildModel(gene, locus, opts={}){
    const spliceRules = Object.assign({ allowSkipEdges:true, requireCanonicalSignals:true, includeU12:false }, opts.spliceRules||{});
    const breakpoints = buildBreakpoints(gene);
    const { segments, segmentIndex, byKey } = buildSegments(gene, breakpoints);
    const transcripts = transcriptPaths(gene, byKey);
    const annotated = buildAnnotatedEdges(transcripts);
    addPlausibleEdges(annotated, { segments, segmentIndex }, locus, spliceRules);
    const edges = Array.from(annotated.values());
    const collapsedPieces = unionExonicPieces(segments);
    const collapsed = { pieces: collapsedPieces, g2c: g2cFactory(collapsedPieces), c2g: c2gFactory(collapsedPieces) };
    const genomic = { start: locus.start, end: locus.end };
    const lanes = assignLanes(segments);
    return { segments, segmentIndex, edges, transcripts, collapsed, genomic, lanes };
  }

  function enumerateChains(model, opts={}){
    const maxChains = opts.maxChains || 200;
    const maxSegments = opts.maxSegments || 12;
    const adj = new Map(); const indeg = new Map(); const outdeg = new Map();
    model.edges.forEach(e => { if (!adj.has(e.from)) adj.set(e.from, []); adj.get(e.from).push(e.to); outdeg.set(e.from, (outdeg.get(e.from)||0)+1); indeg.set(e.to, (indeg.get(e.to)||0)+1); });
    // starts: first segment of each transcript (de-duplicated)
    const startSet = new Set(model.transcripts.map(t => t.segments[0]).filter(Boolean));
    const chains = [];
    function dfs(cur, path){
      if (chains.length >= maxChains) return;
      if (path.length >= maxSegments){ chains.push(path.slice()); return; }
      const nbrs = adj.get(cur) || [];
      if (nbrs.length === 0){ chains.push(path.slice()); return; }
      for (const v of nbrs){ if (path.includes(v)) continue; path.push(v); dfs(v, path); path.pop(); if (chains.length >= maxChains) return; }
    }
    startSet.forEach(s => { if (s){ dfs(s, [s]); } });
    // serialize to exon arrays
    const out = chains.map((segIds, i) => ({ id:`CHAIN_${i+1}`, segments: segIds }));
    return out;
  }

  function assembleChainDNA(model, locus, chain){
    const parts = [];
    (chain.segments||[]).forEach(id => {
      const s = model.segmentIndex.get(id); if (!s) return;
      const dna = orientedSeq(locus, s.start, s.end);
      parts.push(dna);
    });
    return parts.join('');
  }

  function dnaToRNA(s){ return (s||'').toUpperCase().replace(/[^ACGT]/g,'').replace(/T/g,'U'); }

  global.SpliceModel = { buildModel, enumerateChains, assembleChainDNA, dnaToRNA };
})(typeof window !== 'undefined' ? window : globalThis);

