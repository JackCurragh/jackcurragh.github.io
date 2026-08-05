/**
 * SpliceView (SVG + D3)
 * Sashimi-style interactive splicing renderer with collapsed/genomic x-modes.
 */
(function(global){
  const d3 = global.d3;

  function init(container, model, opts={}){
    const width = (opts.width || container.clientWidth || 900);
    const height = (opts.height || 260);
    const margin = { top: 16, right: 16, bottom: 24, left: 16 };
    const laneH = { expanded: 28, dense: 18, collapsed: 12 };
    let viewState = 'expanded';
    let xMode = 'collapsed';
    let intronGapPx = 12;
    // Layout mode:
    //  - 'core' (default): draw each unique block once with flow-guided lane packing
    //  - 'perTranscript': duplicate blocks per transcript lane (teaching/diagnostic)
    //  - 'shared': unique blocks shown once, positioned on the first transcript lane that uses them
    //  - 'graph': true splice graph — one rect per unique exon with independent lane packing
    let layoutMode = 'core';

    container.innerHTML = '';
    const svg = d3.select(container).append('svg').attr('width', width).attr('height', height);
    const defs = svg.append('defs');
    const glow = defs.append('filter').attr('id','glow').attr('height','300%').attr('width','300%').attr('x','-100%').attr('y','-100%');
    glow.append('feGaussianBlur').attr('stdDeviation', '2.0').attr('result','coloredBlur');
    const feMerge = glow.append('feMerge');
    feMerge.append('feMergeNode').attr('in','coloredBlur');
    feMerge.append('feMergeNode').attr('in','SourceGraphic');
    const root = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);
    // Background rect to capture click-away clears. Inserted before content so it doesn't block interactions.
    const bg = root.insert('rect', ':first-child').attr('class','bg')
      .attr('fill','transparent').style('pointer-events','all')
      .on('click', ()=> api._emit('bg-click', {}));
    const content = root.append('g').attr('class','content');
    const axes = root.append('g').attr('class','axes');

    // scales
    const cMax = model.collapsed.pieces.length ? model.collapsed.pieces[model.collapsed.pieces.length-1].cEnd : 1;
    let xCollapsed = d3.scaleLinear().domain([0, cMax]).range([0, width - margin.left - margin.right]);
    let xGenomic = d3.scaleLinear().domain([model.genomic.start, model.genomic.end]).range([0, width - margin.left - margin.right]);
    let x = xCollapsed;
    const pieces = model.collapsed.pieces || [];

    function pieceIndexForGenomic(g){
      // binary search pieces by genomic coord; if intronic, return previous piece index
      let lo=0, hi=pieces.length-1, idx=-1;
      while (lo<=hi){
        const mid=(lo+hi)>>1; const p=pieces[mid];
        if (g < p.start) hi=mid-1; else if (g > p.end) lo=mid+1; else { idx=mid; break; }
      }
      if (idx === -1) idx = Math.max(0, lo-1);
      return Math.max(0, Math.min(pieces.length-1, idx));
    }

    function collapsedWithGaps(g){
      const base = model.collapsed.g2c(g);
      const pi = pieceIndexForGenomic(g);
      return base + (pi * intronGapPx);
    }

    const api = { _handlers: new Map() };
    const zoom = d3.zoom().scaleExtent([0.5, 40]).on('zoom', (ev)=>{
      const t = ev.transform;
      const base = (xMode === 'collapsed') ? xCollapsed : xGenomic;
      x = t.rescaleX(base);
      draw();
      api._emit('view-change', { transform: t });
    });
    svg.call(zoom);
    api.on = (ev, fn) => { api._handlers.set(ev, fn); return api; };
    api._emit = (ev, payload) => { const f = api._handlers.get(ev); if (f) f(payload); };
    api.setMode = ({ xMode:mode, intronGapPx:gap }) => {
      if (mode) {
        xMode = mode;
        // Reset x scale base for the selected mode and clear zoom so users see content immediately
        x = (xMode === 'collapsed') ? xCollapsed : xGenomic;
        try { svg.call(zoom.transform, d3.zoomIdentity); } catch {}
      }
      if (gap!=null) intronGapPx = gap;
      draw();
    };
    api.setViewState = (s) => { viewState = s; draw(); };
    api.highlightChain = (chain) => { highlighted = chain; draw(); };
    api.highlightExon = (exonId) => { highlighted = exonId ? { exonFocus: exonId } : null; draw(); };
    api.setLayoutMode = (mode)=>{ if (mode==='core' || mode==='perTranscript' || mode==='shared' || mode==='graph'){ layoutMode = mode; draw(); } return api; };

    let highlighted = null;

    // Transcript lanes (index by transcript order) for per-transcript layout
    // Sort transcripts by genomic position of their first block to make lanes more readable
    const transcripts = Array.isArray(model.transcripts) ? model.transcripts.slice() : [];
    const blockStart = (bid)=>{ const b=model.blockIndex.get(bid); return b? b.start : Infinity; };
    transcripts.sort((a,b)=> (blockStart(a.blocks?.[0]) - blockStart(b.blocks?.[0])) || String(a.id).localeCompare(String(b.id)) );
    const txIndex = new Map(transcripts.map((t,i)=> [t.id, i]));

    // Exon-aware shared layout: place each unique exon on a transcript lane, not individual blocks
    const exonLaneOf = new Map(); // exonId -> lane index
    const blockToExon = new Map(); // bid -> exonId (primary)
    // Choose the earliest transcript lane that contains the exon
    for (const t of transcripts){
      const lane = txIndex.get(t.id) || 0;
      for (const eid of (t.exons||[])){
        if (!exonLaneOf.has(eid)) exonLaneOf.set(eid, lane);
      }
    }
    // Map blocks to a representative exon id (first covering exon)
    for (const b of (model.blocks||[])){
      const eids = (b.covers && b.covers.exonIds) || [];
      if (eids && eids.length){ blockToExon.set(b.id, eids[0]); }
    }

    // Graph layout lanes: greedy interval packing on unique exon spans (independent of transcript lanes)
    const exonGraphLaneOf = new Map();
    let exonGraphLaneCount = 1;
    (function assignExonGraphLanes(){
      const exons = Array.isArray(model.exons) ? model.exons.slice().sort((a,b)=> a.start - b.start) : [];
      const laneEnds = [];
      for (const ex of exons){
        let lane = -1;
        for (let i=0;i<laneEnds.length;i++){
          if (ex.start >= laneEnds[i]){ lane = i; laneEnds[i] = ex.end; break; }
        }
        if (lane === -1){ lane = laneEnds.length; laneEnds.push(ex.end); }
        exonGraphLaneOf.set(ex.id, lane);
      }
      exonGraphLaneCount = Math.max(1, laneEnds.length||1);
    })();

    function laneOfKey(key){
      // In block-union layout all blocks tile without overlap (-> lane 0).
      // When rendering per transcript, lane is driven by transcript index.
      try { return (model.lanes && model.lanes.laneOf && model.lanes.laneOf.get(key)) || 0; } catch { return 0; }
    }
    function lanesCount(){
      if (layoutMode === 'perTranscript') return Math.max(1, transcripts.length);
      if (layoutMode === 'graph') return Math.max(1, exonGraphLaneCount);
      return model.lanes.laneCount || 1;
    }

    function xForGenomic(g){ return x(g); }
    function xForCollapsed(g){ return x(collapsedWithGaps(g)); }
    function xForSegment(seg){
      const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
      return [ xf(seg.start), xf(seg.end) ];
    }

    function arcPath(a, b){
      const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
      const x1 = xf(a.g), x2 = xf(b.g);
      const y1 = yForLane(laneOfKey(a.id));
      const y2 = yForLane(laneOfKey(b.id));
      // Control point above both lanes; higher if lanes are separated
      const dy = Math.abs(y1 - y2);
      const ctrlY = Math.min(y1, y2) - Math.max(18, 8 + dy * 0.4);
      const mid = (x1 + x2) / 2;
      return `M${x1},${y1} Q ${mid},${ctrlY} ${x2},${y2}`;
    }

    function laneHeight(){ return laneH[viewState] || 18; }
    function rectHeight(){ return Math.max(8, laneHeight() - 6); }
    function yForLane(l){
      const h = laneHeight();
      return 8 + (l * (h + 6)) + 32; // centerline per lane
    }
    function rectTopForLane(l){ return yForLane(l) - rectHeight()/2; }

    function drawAxis(){
      axes.selectAll('*').remove();
      const g = axes.append('g').attr('transform', 'translate(0,12)');
      const axis = d3.axisTop(x).ticks(6);
      g.call(axis);
    }

    function drawIntronsCollapsed(){
      // draw zigzags at boundaries between consecutive pieces (gap visualization)
      const boundaries = [];
      for (let i=0;i<pieces.length-1;i++) boundaries.push({ left: pieces[i], right: pieces[i+1], idx:i });
      const group = content.selectAll('g.intron').data(boundaries, d => d.left.end + ':' + d.right.start);
      const enter = group.enter().append('g').attr('class','intron');
      enter.append('path').attr('stroke','#64748b').attr('fill','none').attr('stroke-width',1.5);
      const all = enter.merge(group);
      all.select('path').attr('d', d => {
        const gx = d.left.end; // boundary genomic coord
        const x1 = x(collapsedWithGaps(gx));
        const y = yForLane(0) - 32; // near top
        const step = Math.max(6, intronGapPx/3);
        // draw across [x1, x1 + intronGapPx]
        return `M${x1},${y} l${step/2},-6 l${step/2},6 l${step/2},-6 l${step/2},6`;
      }).style('display', intronGapPx>0? 'block':'none');
      group.exit().remove();
    }

    function drawSegments(){
      // Draw exon blocks. In perTranscript mode, render a copy of each block on each transcript lane it belongs to.
      const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
      const sel = new Set(nodesList||[]);
      const isExonSelection = Array.isArray(nodesList) && nodesList.length>0 && model.exonIndex && nodesList.every(id=>model.exonIndex.has(id));
      const selBlocksFromExons = new Set();
      if (isExonSelection){
        try {
          nodesList.forEach(eid => { const ex = model.exonIndex.get(eid); (ex?.blocks||[]).forEach(bid => selBlocksFromExons.add(bid)); });
        } catch {}
      }
      const focusExonId = highlighted && highlighted.exonFocus;
      const focusBlocks = new Set();
      if (focusExonId && model.exonIndex && model.exonIndex.has(focusExonId)){
        try { (model.exonIndex.get(focusExonId).blocks||[]).forEach(bid => focusBlocks.add(bid)); } catch {}
      }

      if (layoutMode === 'perTranscript'){
        // Flatten [tx, blockId] pairs
        const rows = [];
        for (const t of transcripts){
          const lane = txIndex.get(t.id) || 0;
          for (const bid of (t.blocks||[])){
            const b = model.blockIndex.get(bid); if (!b) continue;
            rows.push({ key:`${t.id}:${bid}`, txId:t.id, lane, bid, start:b.start, end:b.end });
          }
        }
        const group = content.selectAll('g.span').data(rows, d=> d.key);
        const enter = group.enter().append('g').attr('class','span');
        enter.append('rect').attr('rx',3).attr('ry',3)
          .attr('fill','#0ea5e9').attr('fill-opacity',0.18).attr('stroke','#0ea5e9').attr('stroke-width',1.5)
          .on('click', (_,d)=> api._emit('block-click', { blockId: d.bid, txId: d.txId, exonId: (blockToExon.get(d.bid)||null) }));
        enter.append('title');
        const all = enter.merge(group);
        all.select('rect')
          .attr('x', d=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; return xf(d.start); })
          .attr('y', d=> rectTopForLane(d.lane))
          .attr('width', d=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; const x1=xf(d.start), x2=xf(d.end); return Math.max(2, x2-x1); })
          .attr('height', rectHeight())
          .attr('fill-opacity', d => (sel.has(d.bid) || selBlocksFromExons.has(d.bid) || focusBlocks.has(d.bid)) ? 0.35 : 0.18)
          .attr('stroke', d => (sel.has(d.bid) || selBlocksFromExons.has(d.bid) || focusBlocks.has(d.bid)) ? '#f59e0b' : '#0ea5e9')
          .attr('stroke-width', d => (sel.has(d.bid) || selBlocksFromExons.has(d.bid) || focusBlocks.has(d.bid)) ? 2.5 : 1.5)
          .attr('filter', d => (sel.has(d.bid) || selBlocksFromExons.has(d.bid) || focusBlocks.has(d.bid)) ? 'url(#glow)' : null);
        all.select('title').text(d=> `TX ${d.txId} · ${d.start}-${d.end}`);
        group.exit().remove();
      } else if (layoutMode === 'graph') {
        // True splice graph: one rect per unique exon on independently packed lanes
        const rows = (model.exons||[]).map(ex => ({
          key:`graph-exon:${ex.id}`, eid: ex.id, start: ex.start, end: ex.end, lane: exonGraphLaneOf.get(ex.id) || 0
        }));
        const group = content.selectAll('g.span').data(rows, d=> d.key);
        const enter = group.enter().append('g').attr('class','span');
        enter.append('rect').attr('rx',4).attr('ry',4)
          .attr('fill','#0ea5e9').attr('fill-opacity',0.2).attr('stroke','#0ea5e9').attr('stroke-width',2)
          .on('click', (_,d)=>{
            try {
              const ex = model.exonIndex.get(d.eid);
              const bid = (ex && ex.blocks && ex.blocks[0]) || null;
              api._emit('block-click', { blockId: bid, exonId: d.eid });
            } catch { /* noop */ }
          });
        enter.append('title');
        const all = enter.merge(group);
        all.select('rect')
          .attr('x', d=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; return xf(d.start); })
          .attr('y', d=> rectTopForLane(d.lane))
          .attr('width', d=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; const x1=xf(d.start), x2=xf(d.end); return Math.max(2, x2-x1); })
          .attr('height', rectHeight())
          .attr('fill-opacity', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? 0.35 : 0.2)
          .attr('stroke', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? '#f59e0b' : '#0ea5e9')
          .attr('stroke-width', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? 3 : 2)
          .attr('filter', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? 'url(#glow)' : null);
        all.select('title').text(d=> `EXON ${d.start}-${d.end}`);
        group.exit().remove();
      } else if (layoutMode === 'shared') {
        // Unique exons; draw spans that cover the exon (merge block spans visually) on the chosen transcript lane
        const rows = [];
        for (const ex of (model.exons||[])){
          const lane = exonLaneOf.get(ex.id) || 0;
          // Use exon coordinates directly for the rectangle
          rows.push({ key:`shared-exon:${ex.id}`, eid: ex.id, start: ex.start, end: ex.end, lane });
        }
        const group = content.selectAll('g.span').data(rows, d=> d.key);
        const enter = group.enter().append('g').attr('class','span');
        enter.append('rect').attr('rx',4).attr('ry',4)
          .attr('fill','#0ea5e9').attr('fill-opacity',0.2).attr('stroke','#0ea5e9').attr('stroke-width',2)
          .on('click', (_,d)=> {
            // emit as block-click using the first block of the exon to reuse handlers
            try {
              const ex = model.exonIndex.get(d.eid);
              const bid = (ex && ex.blocks && ex.blocks[0]) || null;
              api._emit('block-click', { blockId: bid, exonId: d.eid });
            } catch { /* noop */ }
          });
        enter.append('title');
        const all = enter.merge(group);
        all.select('rect')
          .attr('x', d=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; return xf(d.start); })
          .attr('y', d=> rectTopForLane(d.lane))
          .attr('width', d=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; const x1=xf(d.start), x2=xf(d.end); return Math.max(2, x2-x1); })
          .attr('height', rectHeight())
          .attr('fill-opacity', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? 0.35 : 0.2)
          .attr('stroke', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? '#f59e0b' : '#0ea5e9')
          .attr('stroke-width', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? 3 : 2)
          .attr('filter', d => (sel.has(d.eid) || (highlighted && highlighted.exonFocus === d.eid)) ? 'url(#glow)' : null);
        all.select('title').text(d=> `EXON ${d.start}-${d.end}`);
        group.exit().remove();

        // Overlay: when focusing a single exon, highlight copies on all transcript lanes that contain it
        const focusId = highlighted && highlighted.exonFocus;
        let focusRows = [];
        if (focusId && model.exonIndex && model.exonIndex.has(focusId)){
          try {
            const ex = model.exonIndex.get(focusId);
            const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic;
            const x1=xf(ex.start), x2=xf(ex.end);
            const txs = Array.isArray(ex.transcriptIds) ? ex.transcriptIds : [];
            focusRows = txs.filter(tid => txIndex.has(tid)).map(tid => ({ key:`focus:${tid}:${focusId}`, lane: txIndex.get(tid), x1, w: Math.max(2, x2-x1) }));
          } catch {}
        }
        const fsel = content.selectAll('rect.shared-focus').data(focusRows, d=> d.key);
        const fenter = fsel.enter().append('rect').attr('class','shared-focus').style('pointer-events','none').attr('rx',4).attr('ry',4)
          .attr('fill','#f59e0b').attr('fill-opacity',0.15).attr('stroke','#f59e0b').attr('stroke-width',2.5).attr('filter','url(#glow)');
        const fall = fenter.merge(fsel);
        fall
          .attr('x', d=> d.x1)
          .attr('y', d=> rectTopForLane(d.lane))
          .attr('width', d=> d.w)
          .attr('height', rectHeight());
        fsel.exit().remove();
      } else {
        // Fallback: union-of-blocks layout (single lane)
        const spans = model.blocks || [];
        const group = content.selectAll('g.span').data(spans, b=> b.id);
        const enter = group.enter().append('g').attr('class','span');
        enter.append('rect').attr('rx',4).attr('ry',4).attr('fill','#0ea5e9').attr('fill-opacity',0.2).attr('stroke','#0ea5e9').attr('stroke-width',2)
          .on('click', (_,b)=> api._emit('block-click', { blockId: b.id, exonId: (blockToExon.get(b.id)||null) }));
        enter.append('title');
        const all = enter.merge(group);
        all.select('rect')
          .attr('x', b=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; return xf(b.start); })
          .attr('y', b=> rectTopForLane(laneOfKey(b.id)))
          .attr('width', b=>{ const xf=(xMode==='collapsed')?xForCollapsed:xForGenomic; const x1=xf(b.start), x2=xf(b.end); return Math.max(2, x2-x1); })
          .attr('height', rectHeight())
          .attr('fill-opacity', b => (sel.has(b.id) || selBlocksFromExons.has(b.id) || focusBlocks.has(b.id)) ? 0.3 : 0.2);
        all.select('title').text(b=> `BLOCK ${b.start}-${b.end}`);
        group.exit().remove();
      }
    }

    function drawJunctions(){
      const baseEdges = (model.edges || []).filter(e => {
        const a=model.blockIndex.get(e.from), b=model.blockIndex.get(e.to);
        if (!a||!b) return false;
        // Half-open semantics: allow touching boundaries
        return (model.strand===1) ? (b.start >= a.end) : (b.end <= a.start);
      });
      if (layoutMode === 'perTranscript'){
        // one arc per (edge × transcript)
        const rows = [];
        for (const e of baseEdges){
          const txs = Array.isArray(e.transcriptIds) && e.transcriptIds.length ? e.transcriptIds : [null];
          for (const txId of txs){
            const lane = (txId && txIndex.has(txId)) ? txIndex.get(txId) : 0;
            rows.push({ key:`${txId||'_'}.${e.from}->${e.to}`, from:e.from, to:e.to, annotated:!!e.annotated, lane, txId });
          }
        }
        const group = content.selectAll('path.junction').data(rows, d=> d.key);
        const enter = group.enter().append('path').attr('class','junction').attr('fill','none')
          .on('click', (_,d)=> api._emit('edge-click', { from: d.from, to: d.to, kind: 'SPLICE' }));
        const all = enter.merge(group);
        const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
        const nodesSet = new Set(nodesList||[]);
        const isExonSelection = Array.isArray(nodesList) && nodesList.length>0 && model.exonIndex && nodesList.every(id=>model.exonIndex.has(id));
        let selectedEdgeKeys = new Set();
        try {
          if (isExonSelection && global.SpliceCore && typeof global.SpliceCore.edgesFromExons==='function'){
            (global.SpliceCore.edgesFromExons(model, nodesList)||[]).forEach(p=> selectedEdgeKeys.add(p.from+'->'+p.to));
          } else if (Array.isArray(nodesList) && global.SpliceCore && typeof global.SpliceCore.edgesFromBlocks==='function'){
            (global.SpliceCore.edgesFromBlocks(model, nodesList)||[]).forEach(p=> selectedEdgeKeys.add(p.from+'->'+p.to));
          }
        } catch {}
        const focusExonId = highlighted && highlighted.exonFocus;
        const focusTx = new Set();
        if (focusExonId && model.exonIndex && model.exonIndex.has(focusExonId)){
          (model.exonIndex.get(focusExonId).transcriptIds||[]).forEach(tid => focusTx.add(tid));
        }
        all.attr('d', d => {
          const a = model.blockIndex.get(d.from), b = model.blockIndex.get(d.to); if (!a||!b) return '';
          const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y = yForLane(d.lane);
          const mid = (x1 + x2) / 2; const ctrlY = y - Math.max(18, 40);
          return `M${x1},${y} Q ${mid},${ctrlY} ${x2},${y}`;
        })
        .attr('stroke', d => (selectedEdgeKeys.has(d.from+'->'+d.to) || focusTx.has(d.txId)) ? '#f59e0b' : (d.annotated ? '#22c55e' : '#94a3b8'))
        .attr('stroke-width', d => (selectedEdgeKeys.has(d.from+'->'+d.to) || focusTx.has(d.txId)) ? 3.5 : (d.annotated ? 2 : 1.5))
        .attr('filter', d => (selectedEdgeKeys.has(d.from+'->'+d.to) || focusTx.has(d.txId)) ? 'url(#glow)' : null)
        .attr('stroke-linecap','round')
        .attr('stroke-dasharray', d => d.annotated ? null : '6,4');
        group.exit().remove();
      } else if (layoutMode === 'graph'){
        // Graph: single arc per exon edge anchored to exon-graph lanes
        const group = content.selectAll('path.junction').data((model.exonEdges||[]), e=> e.from + '->' + e.to);
        const enter = group.enter().append('path').attr('class','junction').attr('fill','none')
          .on('click', (_,e)=> api._emit('edge-click', { from: e.from, to: e.to, kind: 'SPLICE' }));
        const all = enter.merge(group);
        const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
        const selectedEdgeKeys = new Set();
        try {
          if (Array.isArray(nodesList) && model.exonIndex && nodesList.length && nodesList.every(id=>model.exonIndex.has(id)) && global.SpliceCore && typeof global.SpliceCore.edgesFromExons==='function'){
            (global.SpliceCore.edgesFromExons(model, nodesList)||[]).forEach(p=> selectedEdgeKeys.add(p.from+'->'+p.to));
          }
        } catch {}
        const focusExonId = highlighted && highlighted.exonFocus;
        all.attr('d', e => {
          const a = model.exonIndex.get(e.from), b = model.exonIndex.get(e.to); if (!a||!b) return '';
          const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y1 = yForLane(exonGraphLaneOf.get(a.id) || 0);
          const y2 = yForLane(exonGraphLaneOf.get(b.id) || 0);
          const dy = Math.abs(y1 - y2);
          const ctrlY = Math.min(y1, y2) - Math.max(18, 8 + dy * 0.4);
          const mid = (x1 + x2) / 2;
          return `M${x1},${y1} Q ${mid},${ctrlY} ${x2},${y2}`;
        })
        .attr('stroke', e => (selectedEdgeKeys.has(e.from+'->'+e.to) || (focusExonId && (e.from===focusExonId || e.to===focusExonId))) ? '#f59e0b' : (e.annotated ? '#22c55e' : '#94a3b8'))
        .attr('stroke-width', e => (selectedEdgeKeys.has(e.from+'->'+e.to) || (focusExonId && (e.from===focusExonId || e.to===focusExonId))) ? 3.5 : (e.annotated ? 2 : 1.5))
        .attr('filter', e => (selectedEdgeKeys.has(e.from+'->'+e.to) || (focusExonId && (e.from===focusExonId || e.to===focusExonId))) ? 'url(#glow)' : null)
        .attr('stroke-linecap','round')
        .attr('stroke-dasharray', e => e.annotated ? null : '6,4');
        group.exit().remove();
      } else if (layoutMode === 'shared'){
        // exon-level arcs: use model.exonEdges; anchor to exon lanes
        const group = content.selectAll('path.junction').data((model.exonEdges||[]), e=> e.from + '->' + e.to);
        const enter = group.enter().append('path').attr('class','junction').attr('fill','none')
          .on('click', (_,e)=> api._emit('edge-click', { from: e.from, to: e.to, kind: 'SPLICE' }));
        const all = enter.merge(group);
        const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
        const sel = new Set(nodesList||[]);
        const selectedEdgeKeys = new Set();
        try {
          if (Array.isArray(nodesList) && model.exonIndex && nodesList.length && nodesList.every(id=>model.exonIndex.has(id)) && global.SpliceCore && typeof global.SpliceCore.edgesFromExons==='function'){
            (global.SpliceCore.edgesFromExons(model, nodesList)||[]).forEach(p=> selectedEdgeKeys.add(p.from+'->'+p.to));
          }
        } catch {}
        all.attr('d', e => {
          const a = model.exonIndex.get(e.from), b = model.exonIndex.get(e.to); if (!a||!b) return '';
          const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y1 = yForLane(exonLaneOf.get(a.id) || 0);
          const y2 = yForLane(exonLaneOf.get(b.id) || 0);
          const dy = Math.abs(y1 - y2);
          const ctrlY = Math.min(y1, y2) - Math.max(18, 8 + dy * 0.4);
          const mid = (x1 + x2) / 2;
          return `M${x1},${y1} Q ${mid},${ctrlY} ${x2},${y2}`;
        })
        .attr('stroke', e => (selectedEdgeKeys.has(e.from+'->'+e.to)) ? '#f59e0b' : (e.annotated ? '#22c55e' : '#94a3b8'))
        .attr('stroke-width', e => (selectedEdgeKeys.has(e.from+'->'+e.to)) ? 3.5 : (e.annotated ? 2 : 1.5))
        .attr('filter', e => (selectedEdgeKeys.has(e.from+'->'+e.to)) ? 'url(#glow)' : null)
        .attr('stroke-linecap','round')
        .attr('stroke-dasharray', e => e.annotated ? null : '6,4');
        group.exit().remove();
      } else {
        // core: single arc per edge anchored to the lanes of its endpoint blocks
        const group = content.selectAll('path.junction').data(baseEdges, e=> e.from + '->' + e.to);
        const enter = group.enter().append('path').attr('class','junction').attr('fill','none')
          .on('click', (_,e)=> api._emit('edge-click', { from: e.from, to: e.to, kind: 'SPLICE' }));
        const all = enter.merge(group);
        const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
        const sel = new Set(nodesList||[]);
        const selectedEdgeKeys = new Set();
        try {
          if (Array.isArray(nodesList) && model.exonIndex && nodesList.length && nodesList.every(id=>model.exonIndex.has(id)) && global.SpliceCore && typeof global.SpliceCore.edgesFromExons==='function'){
            (global.SpliceCore.edgesFromExons(model, nodesList)||[]).forEach(p=> selectedEdgeKeys.add(p.from+'->'+p.to));
          } else if (Array.isArray(nodesList) && global.SpliceCore && typeof global.SpliceCore.edgesFromBlocks==='function'){
            (global.SpliceCore.edgesFromBlocks(model, nodesList)||[]).forEach(p=> selectedEdgeKeys.add(p.from+'->'+p.to));
          }
        } catch {}
        all.attr('d', e => {
          const a = model.blockIndex.get(e.from), b = model.blockIndex.get(e.to); if (!a||!b) return '';
          const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y1 = yForLane(laneOfKey(a.id));
          const y2 = yForLane(laneOfKey(b.id));
          const dy = Math.abs(y1 - y2);
          const ctrlY = Math.min(y1, y2) - Math.max(18, 8 + dy * 0.4);
          const mid = (x1 + x2) / 2;
          return `M${x1},${y1} Q ${mid},${ctrlY} ${x2},${y2}`;
        })
        .attr('stroke', e => (selectedEdgeKeys.has(e.from+'->'+e.to)) ? '#f59e0b' : (e.annotated ? '#22c55e' : '#94a3b8'))
        .attr('stroke-width', e => (selectedEdgeKeys.has(e.from+'->'+e.to)) ? 3.5 : (e.annotated ? 2 : 1.5))
        .attr('filter', e => (selectedEdgeKeys.has(e.from+'->'+e.to)) ? 'url(#glow)' : null)
        .attr('stroke-linecap','round')
        .attr('stroke-dasharray', e => e.annotated ? null : '6,4');
        group.exit().remove();
      }
    }

    function drawJunctionHighlights(){
      const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
      if (!nodesList || nodesList.length < 2) { content.selectAll('path.junction-hi').remove(); return; }
      // Derive pairs from the model using SpliceCore to ensure we only highlight drawable edges
      let pairs = [];
      try {
        if (global.SpliceCore && typeof global.SpliceCore.edgesFromBlocks === 'function') {
          pairs = global.SpliceCore.edgesFromBlocks(model, nodesList) || [];
        } else {
          // Fallback to naive consecutive pairs
          for (let i=0;i<nodesList.length-1;i++){ pairs.push({ from: nodesList[i], to: nodesList[i+1] }); }
        }
      } catch { for (let i=0;i<nodesList.length-1;i++){ pairs.push({ from: nodesList[i], to: nodesList[i+1] }); } }
      const selected = new Set(pairs.map(p => p.from + '->' + p.to));
      if (layoutMode === 'perTranscript'){
        const rows = [];
        for (const e of (model.edges||[])){
          const key = e.from + '->' + e.to; if (!selected.has(key)) continue;
          const txs = Array.isArray(e.transcriptIds) && e.transcriptIds.length ? e.transcriptIds : [null];
          for (const txId of txs){
            const lane = (txId && txIndex.has(txId)) ? txIndex.get(txId) : 0;
            rows.push({ key:'hi:'+ (txId||'_') + ':' + key, from:e.from, to:e.to, lane });
          }
        }
        const group = content.selectAll('path.junction-hi').data(rows, d=> d.key);
        const enter = group.enter().append('path').attr('class','junction-hi').attr('fill','none').style('pointer-events','none');
        const all = enter.merge(group);
        all.attr('d', d => {
          const a = model.blockIndex.get(d.from), b = model.blockIndex.get(d.to); if (!a||!b) return '';
          const xf = (xMode==='collapsed')?xForCollapsed:xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y = yForLane(d.lane);
          const mid = (x1 + x2) / 2; const ctrlY = y - 60;
          return `M${x1},${y} Q ${mid},${ctrlY} ${x2},${y}`;
        }).attr('stroke', '#f59e0b').attr('stroke-width', 5).attr('stroke-linecap','round').attr('stroke-dasharray', null);
        group.exit().remove();
      } else if (layoutMode === 'shared'){
        // Highlight at exon level
        const rows = [];
        for (const e of (model.exonEdges||[])){
          const key = e.from + '->' + e.to; if (!selected.has(key)) continue;
          rows.push({ key:'hi:'+key, from:e.from, to:e.to });
        }
        const group = content.selectAll('path.junction-hi').data(rows, d=> d.key);
        const enter = group.enter().append('path').attr('class','junction-hi').attr('fill','none').style('pointer-events','none');
        const all = enter.merge(group);
        all.attr('d', d => {
          const a = model.exonIndex.get(d.from), b = model.exonIndex.get(d.to); if (!a||!b) return '';
          const xf = (xMode==='collapsed')?xForCollapsed:xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y1 = yForLane(exonLaneOf.get(a.id) || 0);
          const y2 = yForLane(exonLaneOf.get(b.id) || 0);
          const dy = Math.abs(y1 - y2);
          const ctrlY = Math.min(y1, y2) - Math.max(18, 8 + dy * 0.4);
          const mid = (x1 + x2) / 2;
          return `M${x1},${y1} Q ${mid},${ctrlY} ${x2},${y2}`;
        }).attr('stroke', '#f59e0b').attr('stroke-width', 5).attr('stroke-linecap','round').attr('stroke-dasharray', null);
        group.exit().remove();
      } else {
        const rows = [];
        for (const e of (model.edges||[])){
          const key = e.from + '->' + e.to; if (!selected.has(key)) continue;
          rows.push({ key:'hi:'+key, from:e.from, to:e.to });
        }
        const group = content.selectAll('path.junction-hi').data(rows, d=> d.key);
        const enter = group.enter().append('path').attr('class','junction-hi').attr('fill','none').style('pointer-events','none');
        const all = enter.merge(group);
        all.attr('d', d => {
          const a = model.blockIndex.get(d.from), b = model.blockIndex.get(d.to); if (!a||!b) return '';
          const xf = (xMode==='collapsed')?xForCollapsed:xForGenomic;
          const x1 = xf((model.strand===1) ? a.end : a.start);
          const x2 = xf((model.strand===1) ? b.start : b.end);
          const y1 = yForLane(laneOfKey(a.id));
          const y2 = yForLane(laneOfKey(b.id));
          const dy = Math.abs(y1 - y2);
          const ctrlY = Math.min(y1, y2) - Math.max(18, 8 + dy * 0.4);
          const mid = (x1 + x2) / 2;
          return `M${x1},${y1} Q ${mid},${ctrlY} ${x2},${y2}`;
        }).attr('stroke', '#f59e0b').attr('stroke-width', 5).attr('stroke-linecap','round').attr('stroke-dasharray', null);
        group.exit().remove();
      }
    }

    function drawHighlight(){
      // Highlight selected blocks as translucent bands on their lanes
      content.selectAll('rect.hl').remove();
      const nodesList = Array.isArray(highlighted) ? highlighted : (highlighted && highlighted.nodes) || null;
      if (!nodesList || nodesList.length === 0) return;
      const xf = (xMode === 'collapsed') ? xForCollapsed : xForGenomic;
      // If perTranscript, draw highlight on every transcript lane containing that block
      const rows = [];
      const presentByTx = new Map(); // bid -> Set(txId)
      for (const t of transcripts){ for (const bid of (t.blocks||[])){ if (!presentByTx.has(bid)) presentByTx.set(bid, new Set()); presentByTx.get(bid).add(t.id); } }
      for (const bid of nodesList){
        const b = model.blockIndex.get(bid); if (!b) continue;
        if (layoutMode === 'perTranscript'){
          const set = presentByTx.get(bid) || new Set([null]);
          for (const txId of set){ const lane = txIndex.has(txId) ? txIndex.get(txId) : 0; rows.push({ start:b.start, end:b.end, lane, key:`hl:${txId||'_'}:${bid}` }); }
        } else {
          rows.push({ start:b.start, end:b.end, lane: laneOfKey(b.id), key:`hl:${bid}` });
        }
      }
      const group = content.selectAll('rect.hl').data(rows, d=> d.key);
      const enter = group.enter().append('rect').attr('class','hl').attr('fill','#f59e0b').attr('fill-opacity',0.2).attr('stroke','#f59e0b').attr('stroke-width',2).style('pointer-events','none');
      const all = enter.merge(group);
      all.attr('x', d=> xf(d.start)).attr('y', d=> yForLane(d.lane) - 26).attr('width', d=> Math.max(2, xf(d.end) - xf(d.start))).attr('height', 52);
      group.exit().remove();
    }

    function draw(){
      drawAxis();
      if (xMode === 'collapsed') drawIntronsCollapsed(); else content.selectAll('g.intron').remove();
      // Clear previous elements to avoid stacking artifacts
      content.selectAll('g.span').remove();
      content.selectAll('path.junction').remove();
      content.selectAll('path.junction-hi').remove();
      // Maintain a background rect span that covers the drawing area for click-away clears
      try {
        const w = width - margin.left - margin.right;
        const total = Math.max(1, lanesCount());
        const bottomPad = Math.max(60, Math.round(rectHeight()/2) + 60);
        const h = yForLane(total-1) + bottomPad;
        const finalSvgHeight = Math.max(height, h + margin.top + margin.bottom);
        svg.attr('height', finalSvgHeight);
        bg.attr('x', 0).attr('y', 0).attr('width', w).attr('height', finalSvgHeight);
      } catch {}
      drawSegments();
      drawJunctions();
    }

    draw();
    return api;
  }

  global.SpliceView = { init };
})(typeof window !== 'undefined' ? window : globalThis);
