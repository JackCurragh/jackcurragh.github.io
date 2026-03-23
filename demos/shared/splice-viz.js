/**
 * SpliceViz - Sashimi-style visualization for locus splicing graphs
 *
 * Renders:
 *  - Genomic axis
 *  - Exons as boxes
 *  - Junctions as arcs (thickness ∝ transcript support; dashed if plausible-only)
 *  - Optional highlight of a selected chain (exon route)
 *
 * Usage:
 *   SpliceViz.draw(ctx, canvas, eventModel, locus, { highlightChain })
 */
(function(global){
  function ntToPxFactory(locus, canvas, margin){
    const width = canvas.width - margin * 2;
    const L = Math.max(1, locus.end - locus.start);
    return (g) => margin + ((g - locus.start) / L) * width;
  }

  function resolveNodeKey(node){
    return `${node.id || (node.start+'-'+node.end)}:${node.strand}`;
  }

  function parseEdgeKey(edgeKey){
    // edgeKey like "ENSE...:1->ENSE...:1" or "39520-39600:1->..."
    const [from, to] = String(edgeKey).split('->');
    const parse = (k)=>{
      const part = String(k).split(':')[0];
      if (/^\d+-\d+$/.test(part)){
        const [s,e] = part.split('-').map(x=>parseInt(x,10));
        return { start:s, end:e };
      }
      return { id: part };
    };
    return { from: parse(from), to: parse(to) };
  }

  function findNodeByRef(nodes, ref){
    if (ref.id){
      const n = nodes.find(n=> String(n.id) === String(ref.id));
      if (n) return n;
    }
    if (typeof ref.start === 'number' && typeof ref.end === 'number'){
      // exact match by coords
      const n = nodes.find(n=> n.start === ref.start && n.end === ref.end);
      if (n) return n;
      // nearest by overlap
      let best=null, bestOv=-1;
      for (const n2 of nodes){
        const ov = Math.min(n2.end, ref.end) - Math.max(n2.start, ref.start);
        if (ov>bestOv){ best=n2; bestOv=ov; }
      }
      return best;
    }
    return null;
  }

  function buildEdgeIndex(eventModel){
    const nodes = eventModel.nodes || [];
    const edges = eventModel.edges || [];
    const map = new Map();
    let maxSupport = 1;
    for (const e of edges){
      const key = `${e.from}->${e.to}`;
      const ref = parseEdgeKey(key);
      const n1 = findNodeByRef(nodes, ref.from);
      const n2 = findNodeByRef(nodes, ref.to);
      if (!n1 || !n2) continue;
      const support = Array.isArray(e.transcripts) ? e.transcripts.length : 0;
      maxSupport = Math.max(maxSupport, support);
      map.set(key, {
        key,
        from: n1,
        to: n2,
        support,
        plausible: !Array.isArray(e.transcripts) || e.transcripts.length === 0
      });
    }
    return { map, maxSupport };
  }

  function drawAxis(ctx, canvas, margin, y){
    ctx.strokeStyle = '#334155';
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(margin, y); ctx.lineTo(canvas.width - margin, y); ctx.stroke();
  }

  function drawExons(ctx, toX, y, nodes, opts){
    for (const n of nodes){
      const x1 = toX(n.start), x2 = toX(n.end);
      const w = Math.max(2, x2 - x1);
      ctx.fillStyle = opts.exonFill || 'rgba(14,165,233,0.25)';
      ctx.fillRect(x1, y-22, w, 44);
      ctx.strokeStyle = opts.exonStroke || '#0ea5e9';
      ctx.lineWidth = 2;
      ctx.strokeRect(x1, y-22, w, 44);
    }
  }

  function drawArc(ctx, x1, x2, y, h, width, color, dashed){
    const mid = (x1 + x2) / 2;
    ctx.save();
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    if (dashed) ctx.setLineDash([6,4]);
    ctx.beginPath();
    ctx.moveTo(x1, y-22);
    ctx.quadraticCurveTo(mid, y-22-h, x2, y-22);
    ctx.stroke();
    ctx.restore();
  }

  function drawJunctions(ctx, toX, y, edgeIndex){
    // Draw in order of span for nicer layering
    const list = Array.from(edgeIndex.map.values()).sort((a,b)=>{
      const sa = Math.abs(a.to.start - a.from.end);
      const sb = Math.abs(b.to.start - b.from.end);
      return sa - sb;
    });
    const maxSpan = list.reduce((m,e)=> Math.max(m, Math.abs(e.to.start - e.from.end)), 1);
    for (const e of list){
      const x1 = toX(e.from.end);
      const x2 = toX(e.to.start);
      const span = Math.abs(e.to.start - e.from.end);
      const h = Math.max(18, 60 * (span / (maxSpan||1)));
      const w = e.plausible ? 1.5 : 1 + 6 * (e.support / (edgeIndex.maxSupport||1));
      const color = e.plausible ? '#94a3b8' : '#22c55e';
      drawArc(ctx, x1, x2, y, h, w, color, e.plausible);
    }
  }

  function drawChainHighlight(ctx, toX, y, chain){
    if (!chain || !Array.isArray(chain.exons) || chain.exons.length === 0) return;
    // highlight exons
    ctx.save();
    ctx.globalAlpha = 0.35;
    ctx.fillStyle = '#f59e0b';
    for (const ex of chain.exons){
      const x1 = toX(ex.start), x2 = toX(ex.end);
      ctx.fillRect(x1, y-26, Math.max(2, x2-x1), 52);
    }
    ctx.restore();
    // highlight arcs along chain
    ctx.save();
    ctx.strokeStyle = '#f59e0b';
    ctx.lineWidth = 4;
    for (let i=0;i<chain.exons.length-1;i++){
      const a = chain.exons[i];
      const b = chain.exons[i+1];
      const x1 = toX(a.end), x2 = toX(b.start);
      const mid = (x1+x2)/2; const h=70;
      ctx.beginPath(); ctx.moveTo(x1, y-26); ctx.quadraticCurveTo(mid, y-26-h, x2, y-26); ctx.stroke();
    }
    ctx.restore();
  }

  function draw(ctx, canvas, eventModel, locus, opts={}){
    const y = Math.floor(canvas.height/2);
    const margin = 40;
    const toX = ntToPxFactory(locus, canvas, margin);

    ctx.clearRect(0,0,canvas.width,canvas.height);
    drawAxis(ctx, canvas, margin, y);
    drawExons(ctx, toX, y, eventModel.nodes||[], opts);
    const edgeIndex = buildEdgeIndex(eventModel);
    drawJunctions(ctx, toX, y, edgeIndex);
    if (opts.highlightChain) drawChainHighlight(ctx, toX, y, opts.highlightChain);
  }

  global.SpliceViz = { draw };
})(typeof window !== 'undefined' ? window : globalThis);

