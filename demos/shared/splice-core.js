/**
 * SpliceCore - Exon-block splicing core model
 *
 * Core graph: unique exon blocks (nodes) + splice junctions (edges).
 * Transcripts are paths through the block graph.
 */
(function(global){
  function revComp(dna){ const c={A:'T',T:'A',C:'G',G:'C',N:'N'}; return (dna||'').split('').reverse().map(b=>c[b]||'N').join(''); }
  function toRNA(dna){ return (dna||'').toUpperCase().replace(/[^ACGT]/g,'').replace(/T/g,'U'); }
  // Half-open genomic slicing: [gStart, gEndExcl)
  // Returns DNA oriented to strand (reverse-complemented if -1).
  function orientedSlice(locus, gStart, gEndExcl, strand){
    // normalize/clamp to locus; treat end as exclusive
    const s = Math.max(locus.start, Math.min(gStart, gEndExcl));
    const e = Math.min(locus.end + 1, Math.max(gStart, gEndExcl)); // allow exclusive end at locus.end+1
    if (e <= s) return '';
    const so = s - locus.start, eo = e - locus.start; // substring end is exclusive already
    const slice = (locus.dna||'').substring(so, eo);
    return strand === -1 ? revComp(slice) : slice;
  }

  function strandOf(obj, fallback){ return (typeof obj.strand === 'number' && (obj.strand===1 || obj.strand===-1)) ? obj.strand : (fallback||1); }

  function buildBreakpoints(gene){
    const bp = new Set();
    (gene.Transcript||[]).forEach(tr => {
      (tr.Exon||[]).forEach(ex => { bp.add(ex.start); bp.add(ex.end); });
    });
    const arr = Array.from(bp).sort((a,b)=>a-b);
    return arr;
  }

  function buildBlocks(gene, breakpoints){
    // Any span [bp[i], bp[i+1]) covered by ≥1 exon becomes a block
    const blocks = [];
    function covered(s,e){
      for (const tr of (gene.Transcript||[])){
        for (const ex of (tr.Exon||[])){
          const xs = Math.min(ex.start, ex.end), xe = Math.max(ex.start, ex.end);
          if (s >= xs && e <= xe) return true;
        }
      }
      return false;
    }
    for (let i=0;i<breakpoints.length-1;i++){
      const s = breakpoints[i], e = breakpoints[i+1];
      if (s === e) continue;
      if (covered(s,e)) {
        blocks.push({ id:`B${blocks.length+1}`, start:s, end:e, strand: strandOf(gene, 1), covers:{exonIds:[], transcriptIds:[]} });
      }
    }
    const blockIndex = new Map(blocks.map(b => [b.id, b]));
    return { blocks, blockIndex };
  }

  function mapExonToBlocks(blocks, exon){
    const xs = Math.min(exon.start, exon.end), xe = Math.max(exon.start, exon.end);
    const inside = [];
    for (const b of blocks){ if (b.start >= xs && b.end <= xe) inside.push(b); }
    return inside;
  }

  function assignLanes(blocks){
    // Interval graph greedy packing: minimize number of lanes such that
    // blocks with overlapping x-intervals do not share a lane.
    // Since our blocks are split at all shared boundaries, they rarely
    // overlap; this will collapse to as few lanes as possible.
    const spans = blocks.slice().sort((a,b)=> a.start - b.start);
    const laneEnds = [];
    const laneOf = new Map();
    for (const b of spans){
      let placed = false;
      for (let i=0;i<laneEnds.length;i++){
        if (b.start >= laneEnds[i]){ laneEnds[i] = b.end; laneOf.set(b.id, i); placed = true; break; }
      }
      if (!placed){ laneEnds.push(b.end); laneOf.set(b.id, laneEnds.length - 1); }
    }
    return { laneOf, laneCount: laneEnds.length || 1 };
  }

  function buildCoreModel(gene, locus, opts={}){
    const chrom = gene.seq_region_name || gene.seq_region || gene.seq_region_id || 'chr';
    const gStrand = strandOf(gene, 1);
    const breakpoints = buildBreakpoints(gene);
    const { blocks, blockIndex } = buildBlocks(gene, breakpoints);
    const edgesMap = new Map(); // key from->to
    const transcripts = [];
    // Exon aggregates
    const exonMap = new Map(); // exonId -> {id,start,end,strand, transcriptIds:Set, blocks:[]}
    const exonEdgeMap = new Map(); // key exonIdA->exonIdB -> {from,to, annotated:true, transcriptIds:Set}

    (gene.Transcript||[]).forEach(tr => {
      const tStrand = strandOf(tr, gStrand);
      const exons = (tr.Exon||[]).slice().sort((a,b)=> tStrand===1 ? (a.start-b.start) : (b.start-a.start));
      const tBlocks = [];
      const tEdges = [];
      const tExons = [];
      for (let i=0;i<exons.length;i++){
        const ex = exons[i];
        // Register exon aggregate
        if (ex && ex.id){
          if (!exonMap.has(ex.id)){
            // Compute blocks inside exon once (gene orientation independent)
            const insideBlocks = mapExonToBlocks(blocks, ex).sort((a,b)=> a.start - b.start).map(b=> b.id);
            exonMap.set(ex.id, { id: ex.id, start: Math.min(ex.start, ex.end), end: Math.max(ex.start, ex.end), strand: tStrand, transcriptIds: new Set([tr.id]), blocks: insideBlocks });
          } else {
            exonMap.get(ex.id).transcriptIds.add(tr.id);
          }
          tExons.push(ex.id);
        }
        const inside = mapExonToBlocks(blocks, ex);
        const list = (tStrand===1) ? inside.sort((a,b)=>a.start-b.start) : inside.sort((a,b)=>b.start-a.start);
        list.forEach(b => { if (!tBlocks.includes(b.id)) tBlocks.push(b.id); b.covers.transcriptIds.push(tr.id); b.covers.exonIds.push(ex.id); });
        if (i < exons.length-1){
          // splice from last block of this exon to first block of next exon in orientation
          const nextInside = mapExonToBlocks(blocks, exons[i+1]);
          const nextList = (tStrand===1) ? nextInside.sort((a,b)=>a.start-b.start) : nextInside.sort((a,b)=>b.start-a.start);
          const fromBlock = (tStrand===1) ? list[list.length-1] : list[list.length-1];
          const toBlock = (tStrand===1) ? nextList[0] : nextList[0];
          if (fromBlock && toBlock){
            const key = `${fromBlock.id}->${toBlock.id}`;
            if (!edgesMap.has(key)){
              const donorG = (tStrand===1) ? fromBlock.end : fromBlock.start;
              const acceptorG = (tStrand===1) ? toBlock.start : toBlock.end;
              edgesMap.set(key, { from: fromBlock.id, to: toBlock.id, donorG, acceptorG, transcriptIds:[tr.id], annotated:true });
            } else {
              const e = edgesMap.get(key); if (!e.transcriptIds.includes(tr.id)) e.transcriptIds.push(tr.id);
            }
            tEdges.push({ from: fromBlock.id, to: toBlock.id });
            // exon-level edge
            if (ex && ex.id && exons[i+1] && exons[i+1].id){
              const ek = `${ex.id}->${exons[i+1].id}`;
              if (!exonEdgeMap.has(ek)) exonEdgeMap.set(ek, { from: ex.id, to: exons[i+1].id, annotated:true, transcriptIds: new Set() });
              exonEdgeMap.get(ek).transcriptIds.add(tr.id);
            }
          }
        }
      }
      transcripts.push({ id: tr.id, name: tr.display_name || tr.id, strand: tStrand, blocks: tBlocks, exons: tExons, edges: tEdges });
    });

    // Adjacency
    const edges = Array.from(edgesMap.values());
    const adj = new Map(); edges.forEach(e => { if (!adj.has(e.from)) adj.set(e.from, []); adj.get(e.from).push(e.to); });

    // Collapsed pieces: union of block intervals
    const merged = [];
    blocks.slice().sort((a,b)=>a.start-b.start).forEach(b => {
      if (!merged.length || b.start > merged[merged.length-1][1]) merged.push([b.start, b.end]);
      else merged[merged.length-1][1] = Math.max(merged[merged.length-1][1], b.end);
    });
    let c=0; const pieces=[];
    for (const [s,e] of merged){ const len = Math.max(0, e-s); pieces.push({ start:s, end:e, cStart:c, cEnd:c+len }); c += len; }
    function g2c(g){ for (const p of pieces){ if (g>=p.start && g<=p.end) return p.cStart + (g-p.start); } let bestC=0, bestD=Infinity; for(const p of pieces){ if (g<p.start){ const d=p.start-g; if(d<bestD){bestD=d; bestC=p.cStart;} } else if (g>p.end){ const d=g-p.end; if(d<bestD){bestD=d; bestC=p.cEnd;} } } return bestC; }
    function c2g(c){ for (const p of pieces){ if (c>=p.cStart && c<=p.cEnd) return p.start + (c-p.cStart); } return (pieces[0]?.start)||0; }

    const lanes = assignLanes(blocks);
    // Serialize exons
    const exons = Array.from(exonMap.values()).map(e => ({ id:e.id, start:e.start, end:e.end, strand:gStrand, transcriptIds: Array.from(e.transcriptIds), blocks: e.blocks }));
    const exonIndex = new Map(exons.map(e => [e.id, e]));
    const exonEdges = Array.from(exonEdgeMap.values()).map(e => ({ from:e.from, to:e.to, annotated:true, transcriptIds: Array.from(e.transcriptIds) }));
    const genomic = { start: locus.start, end: locus.end };
    return { chrom: chrom, strand: gStrand, genomic, blocks, blockIndex, edges, adj, transcripts, exons, exonIndex, exonEdges, collapsed:{ pieces, g2c, c2g }, lanes, meta: { coordSemantics: 'half-open' } };
  }

  function enumerateChains(model, opts={}){
    const maxChains = opts.maxChains || 200;
    const maxSteps = opts.maxSteps || 12;
    const starts = new Set();
    (model.transcripts||[]).forEach(t => { if (t.blocks && t.blocks.length) starts.add(t.blocks[0]); });
    const chains=[];
    function dfs(nid, path, steps){
      if (chains.length>=maxChains) return;
      if (steps>=maxSteps) { chains.push({ id:`CHAIN_${chains.length+1}`, blocks: path.slice(), edges: edgesFromBlocks(model, path) }); return; }
      const nbrs = model.adj.get(nid) || [];
      if (!nbrs.length){ chains.push({ id:`CHAIN_${chains.length+1}`, blocks: path.slice(), edges: edgesFromBlocks(model, path) }); return; }
      for (const v of nbrs){ if (path.includes(v)) continue; path.push(v); dfs(v, path, steps+1); path.pop(); if (chains.length>=maxChains) return; }
    }
    Array.from(starts).forEach(s => dfs(s, [s], 0));
    return chains;
  }

  function edgesFromBlocks(model, blocks){
    const out=[];
    for (let i=0;i<blocks.length-1;i++){
      const a=blocks[i], b=blocks[i+1];
      const has = (model.edges||[]).find(e=>e.from===a && e.to===b);
      if (has) out.push({from:a,to:b});
    }
    return out;
  }

  function edgesFromExons(model, exons){
    const out=[];
    for (let i=0;i<exons.length-1;i++){
      const a=exons[i], b=exons[i+1];
      const has = (model.exonEdges||[]).find(e=>e.from===a && e.to===b);
      if (has) out.push({from:a,to:b});
    }
    return out;
  }

  function blocksToExons(model, blocks){
    const seen=[]; let last=null;
    for (const bid of blocks){
      const b = model.blockIndex.get(bid); if (!b) continue;
      const eids = (b.covers && b.covers.exonIds) || [];
      const eid = eids && eids.length ? eids[0] : null; // heuristic: first annotated exon
      if (eid && eid !== last){ seen.push(eid); last = eid; }
    }
    return seen;
  }

  function assembleChainDNA(model, locus, chain){
    if (!chain || !Array.isArray(chain.blocks) || chain.blocks.length===0) return '';
    const strand = model.strand || 1;
    const parts=[];
    for (const bid of chain.blocks){
      const b = model.blockIndex.get(bid); if (!b) continue;
      // Half-open [start, end): pass end as exclusive boundary
      const dna = orientedSlice(locus, b.start, b.end, strand);
      parts.push(dna);
    }
    return parts.join('');
  }

  // Public API: expose exon helpers for congruent highlighting
  global.SpliceCore = { buildCoreModel, enumerateChains, assembleChainDNA, toRNA, edgesFromBlocks, edgesFromExons, blocksToExons };
})(typeof window !== 'undefined' ? window : globalThis);
