/**
 * SpliceEvents - Donor/Acceptor event graph and assembly
 *
 * buildEventModel(gene, locus, opts) -> model
 * enumerateChains(model, opts) -> chains
 * assembleChainDNA(model, locus, chain) -> DNA
 * toRNA(dna)
 */
(function(global){
  function revComp(dna){ const c={A:'T',T:'A',C:'G',G:'C',N:'N'}; return (dna||'').split('').reverse().map(b=>c[b]||'N').join(''); }
  function dnaToRNA(s){ return (s||'').toUpperCase().replace(/[^ACGT]/g,'').replace(/T/g,'U'); }
  function orientedSlice(locus, g1, g2){
    const s = Math.max(locus.start, Math.min(g1,g2));
    const e = Math.min(locus.end, Math.max(g1,g2));
    if (e < s) return '';
    const so = s - locus.start, eo = e - locus.start;
    const slice = (locus.dna||'').substring(so, eo+1);
    return locus.strand === -1 ? revComp(slice) : slice;
  }
  function donorIsGT(locus, donorG){
    if (locus.strand === 1){ const s = orientedSlice(locus, donorG+1, donorG+2); return s === 'GT'; }
    const s = orientedSlice(locus, donorG-2, donorG-1); return s === 'GT';
  }
  function acceptorIsAG(locus, acceptorG){
    if (locus.strand === 1){ const s = orientedSlice(locus, acceptorG-2, acceptorG-1); return s === 'AG'; }
    const s = orientedSlice(locus, acceptorG+1, acceptorG+2); return s === 'AG';
  }
  function strandOf(obj, fallback){ return (typeof obj.strand === 'number' && (obj.strand===1 || obj.strand===-1)) ? obj.strand : (fallback||1); }
  function sortExons(tr, geneStrand){
    const s = strandOf(tr, geneStrand);
    const list = (tr.Exon||[]).slice().sort((a,b)=> s===1 ? (a.start-b.start) : (b.start-a.start));
    return { strand:s, exons:list };
  }
  function nodeId(type, chrom, g, strand){ return `${type}|${chrom}|${g}|${strand}`; }
  function buildEventModel(gene, locus, opts={}){
    const chrom = gene.seq_region_name || gene.seq_region || gene.seq_region_id || 'chr';
    const gStrand = strandOf(gene, 1);
    const nodes = new Map();
    const edges = new Map();
    const nodeIndex = new Map();
    const transcripts = [];
    const cdsMarkers = [];

    function ensureNode(type, g, strand, exonId){
      const id = nodeId(type, chrom, g, strand);
      if (!nodes.has(id)) nodes.set(id, { id, type, g, chrom, strand, exonId });
      return nodes.get(id);
    }
    function addEdge(fromId, toId, kind, trId){
      const k = `${fromId}->${toId}`;
      if (!edges.has(k)) edges.set(k, { from: fromId, to: toId, kind, transcriptIds: [], annotated: true });
      const e = edges.get(k);
      if (trId && !e.transcriptIds.includes(trId)) e.transcriptIds.push(trId);
    }

    // Build nodes/edges per transcript
    (gene.Transcript||[]).forEach(tr => {
      const { strand, exons } = sortExons(tr, gStrand);
      if (!exons.length) return;
      const chain = [];
      const edgeSeq = [];
      exons.forEach(ex => {
        const start = ex.start, end = ex.end;
        const A = (strand===1) ? start : end;
        const D = (strand===1) ? end : start;
        const aN = ensureNode('ACCEPTOR', A, strand, ex.id);
        const dN = ensureNode('DONOR', D, strand, ex.id);
        addEdge(aN.id, dN.id, 'EXON_SPAN', tr.id);
        const eK = `${aN.id}->${dN.id}`;
        edges.get(eK).lengthNt = Math.abs(D - A) + 1;
        // Always push A then D so chain order is A1,D1,A2,D2,...
        chain.push(aN.id);
        chain.push(dN.id);
        edgeSeq.push({ from:aN.id, to:dN.id, kind:'EXON_SPAN' });
      });
      // Splice edges between consecutive exons
      for (let i=0;i<exons.length-1;i++){
        const cur = exons[i], nxt = exons[i+1];
        const D = (strand===1) ? cur.end : cur.start;
        const A = (strand===1) ? nxt.start : nxt.end;
        const dId = nodeId('DONOR', chrom, D, strand);
        const aId = nodeId('ACCEPTOR', chrom, A, strand);
        addEdge(dId, aId, 'SPLICE', tr.id);
        // Do not push aId; it's already in chain from the exon loop at position 2*(i+1)
        edgeSeq.push({ from:dId, to:aId, kind:'SPLICE' });
      }
      // Optional boundary nodes
      if (opts.boundaryNodes!==false){
        const tstart = { id: nodeId('TSTART', chrom, chain.length? nodes.get(chain[0]).g : (gene.start||locus.start), strand), type:'TSTART', g: (chain.length? nodes.get(chain[0]).g : (gene.start||locus.start)), chrom, strand };
        const tend = { id: nodeId('TEND', chrom, chain.length? nodes.get(chain[chain.length-1]).g : (gene.end||locus.end), strand), type:'TEND', g: (chain.length? nodes.get(chain[chain.length-1]).g : (gene.end||locus.end)), chrom, strand };
        nodes.set(tstart.id, tstart); nodes.set(tend.id, tend);
        addEdge(tstart.id, chain[0], 'BOUNDARY', tr.id);
        addEdge(chain[chain.length-1], tend.id, 'BOUNDARY', tr.id);
        chain.unshift(tstart.id); chain.push(tend.id);
      }
      // CDS markers
      if (opts.includeCDSMarkers!==false && tr.Translation){
        cdsMarkers.push({ transcriptId: tr.id, startG: tr.Translation.start, stopG: tr.Translation.end });
      }
      transcripts.push({ id: tr.id, name: tr.display_name || tr.id, strand, chain, edges: edgeSeq });
    });

    // Serialize nodes
    const nodeArr = Array.from(nodes.values());
    nodeArr.forEach(n => nodeIndex.set(n.id, n));

    // Collapsed pieces (exonic union)
    const spans = Array.from(edges.values()).filter(e => e.kind==='EXON_SPAN').map(e => {
      const a = nodeIndex.get(e.from), d = nodeIndex.get(e.to); return [Math.min(a.g,d.g), Math.max(a.g,d.g)];
    }).sort((x,y)=> x[0]-y[0]);
    const merged = [];
    for (const [s,e] of spans){
      if (!merged.length || s > merged[merged.length-1][1]) merged.push([s,e]);
      else merged[merged.length-1][1] = Math.max(merged[merged.length-1][1], e);
    }
    let c=0; const pieces=[];
    for (const [s,e] of merged){ const len = Math.max(0, e-s); pieces.push({ start:s, end:e, cStart:c, cEnd:c+len }); c += len; }
    function g2c(g){
      for (const p of pieces){ if (g>=p.start && g<=p.end) return p.cStart + (g - p.start); }
      // intron: clamp to nearest boundary
      let best=null, bestD=Infinity, bestC=0; for (const p of pieces){ if (g < p.start){ const d=p.start-g; if (d<bestD){bestD=d; bestC=p.cStart;} } else if (g>p.end){ const d=g-p.end; if (d<bestD){bestD=d; bestC=p.cEnd;} } }
      return bestC;
    }
    function c2g(c){ for (const p of pieces){ if (c>=p.cStart && c<=p.cEnd) return p.start + (c - p.cStart); } return (pieces[0]?.start)||0; }

    // Adjacency
    const adj = new Map();
    Array.from(edges.values()).forEach(e => { if (!adj.has(e.from)) adj.set(e.from, []); adj.get(e.from).push(e.to); });

    // Lanes: assign based on EXON_SPAN intervals (greedy interval coloring)
    function assignLanes(){
      const spans = Array.from(edges.values()).filter(e => e.kind==='EXON_SPAN').map(e => {
        const a = nodeIndex.get(e.from), d = nodeIndex.get(e.to);
        return { from:e.from, to:e.to, start: Math.min(a.g,d.g), end: Math.max(a.g,d.g) };
      }).sort((x,y)=> x.start - y.start);
      const laneEnds = []; const laneOf = new Map();
      spans.forEach(s => {
        let placed = false;
        for (let i=0;i<laneEnds.length;i++){
          if (s.start > laneEnds[i]){ laneEnds[i] = s.end; laneOf.set(s.from, i); laneOf.set(s.to, i); placed=true; break; }
        }
        if (!placed){ laneEnds.push(s.end); const idx = laneEnds.length-1; laneOf.set(s.from, idx); laneOf.set(s.to, idx); }
      });
      return { laneOf, laneCount: laneEnds.length || 1 };
    }

    const lanes = assignLanes();
    const genomic = { start: locus.start, end: locus.end };
    return { chrom, strand: gStrand, genomic, nodes: nodeArr, nodeIndex, edges: Array.from(edges.values()), adj, transcripts, cdsMarkers, collapsed:{ pieces, g2c, c2g }, lanes };
  }

  function isDownstream(aNode, bNode){
    if (!aNode || !bNode) return false;
    return aNode.strand === 1 ? (bNode.g > aNode.g) : (bNode.g < aNode.g);
  }
  function enumerateChains(model, opts={}){
    const maxChains = opts.maxChains || 200;
    const maxSteps = opts.maxSteps || 12;
    const allowPlausible = (opts.allowPlausible!==false);
    const rules = Object.assign({ requireGTAG:true, includeU12:false }, opts.plausibleRules||{});
    const starts = new Set();
    (model.transcripts||[]).forEach(t => { if (t.chain && t.chain.length){ // pick first ACCEPTOR or first node
      const first = t.chain.find(nid => (model.nodeIndex.get(nid)||{}).type==='ACCEPTOR') || t.chain[0]; starts.add(first);
    }});

    const chains=[];
    const visiting=new Set();
    function neighbors(nid){
      const list = (model.adj.get(nid)||[]).slice();
      if (allowPlausible){
        const from = model.nodeIndex.get(nid);
        if (from && from.type==='DONOR'){
          // plausible to any downstream ACCEPTOR with GT-AG
          for (const node of model.nodes){
            if (node.type!=='ACCEPTOR') continue;
            if (!isDownstream(from, node)) continue;
            if (rules.requireGTAG && !(donorIsGT({start:model.genomic.start,end:model.genomic.end,strand:model.strand,dna:global.__LOCUS_DNA__||''}, from.g) && acceptorIsAG({start:model.genomic.start,end:model.genomic.end,strand:model.strand,dna:global.__LOCUS_DNA__||''}, node.g))) continue;
            if (!list.includes(node.id)) list.push(node.id);
          }
        }
      }
      return list;
    }
    function dfs(nid, path, steps){
      if (chains.length>=maxChains) return;
      if (steps>=maxSteps) { chains.push({ id:`CHAIN_${chains.length+1}`, nodes: path.slice(), edges: [] }); return; }
      const nbrs = neighbors(nid);
      if (!nbrs.length) { chains.push({ id:`CHAIN_${chains.length+1}`, nodes: path.slice(), edges: [] }); return; }
      for (const v of nbrs){ if (path.includes(v)) continue; path.push(v); dfs(v, path, steps+1); path.pop(); if (chains.length>=maxChains) return; }
    }
    Array.from(starts).forEach(s => dfs(s, [s], 0));
    return chains;
  }

  function assembleChainDNA(model, locus, chain){
    if (!chain || !Array.isArray(chain.nodes) || chain.nodes.length<2) return '';
    const parts=[];
    for (let i=0;i<chain.nodes.length-1;i++){
      const a = chain.nodes[i], b = chain.nodes[i+1];
      const k = `${a}->${b}`;
      const e = model.edges.find(e => e.from===a && e.to===b && e.kind==='EXON_SPAN');
      if (!e) continue;
      const an = model.nodeIndex.get(a), dn = model.nodeIndex.get(b);
      const g1 = an.g, g2 = dn.g;
      const dna = orientedSlice(locus, g1, g2);
      parts.push(dna);
    }
    return parts.join('');
  }

  global.SpliceEvents = { buildEventModel, enumerateChains, assembleChainDNA, toRNA: dnaToRNA };
})(typeof window !== 'undefined' ? window : globalThis);
