/**
 * RDG Figure Composer
 *
 * Creates a multi-panel, publication-style figure (A–D) that includes:
 *  A) RDG overview for the focused construct
 *  B) Western blot lanes (saved constructs + current)
 *  C) Construct map for the focused construct
 *  D) Key metrics table for the focused construct
 *
 * Exports a single high-resolution PNG suitable for papers or slides.
 */
(function(root){
  function px(v){ return Math.round(v); }

  // Minimal gel drawing for figure (self-contained)
  function drawGelPanel(ctx, x, y, w, h, lanes){
    // Background
    ctx.save();
    ctx.fillStyle = '#111';
    ctx.fillRect(x, y, w, h);

    const gelTop = y + 50;
    const gelBottom = y + h - 80;
    const gelHeight = gelBottom - gelTop;
    const laneWidth = Math.max(50, Math.min(90, Math.floor(w / Math.max(3, lanes.length+2))));
    const laneGap = Math.max(20, Math.floor((w - laneWidth*lanes.length) / (lanes.length+1)));
    const startX = x + laneGap;
    const minMW = 10, maxMW = 200;

    // MW markers
    ctx.fillStyle = '#bbb';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'right';
    [200,150,100,75,50,37,25,20,15,10].forEach(mw => {
      let frac = (Math.log(mw) - Math.log(minMW)) / (Math.log(maxMW) - Math.log(minMW));
      frac = Math.max(0, Math.min(1, frac));
      const yy = gelTop + gelHeight * (1 - frac);
      ctx.fillText(String(mw), x + 30, yy + 3);
      ctx.strokeStyle = 'rgba(255,255,255,0.08)';
      ctx.beginPath(); ctx.moveTo(x+40, yy); ctx.lineTo(x+w-10, yy); ctx.stroke();
    });

    // Helper to map MW to y
    const yForMW = (mw) => {
      const frac = (Math.log(mw) - Math.log(minMW)) / (Math.log(maxMW) - Math.log(minMW));
      return gelTop + gelHeight * (1 - Math.max(0, Math.min(1, frac)));
    };

    lanes.forEach((lane, i) => {
      const cx = startX + i * (laneWidth + laneGap) + laneWidth/2;
      // Lane background
      ctx.fillStyle = '#2d2d2d';
      ctx.fillRect(cx - laneWidth/2, gelTop, laneWidth, gelHeight);

      const prods = Array.isArray(lane.products) ? lane.products : [];
      const maxAb = Math.max(0.0001, ...prods.map(p => (typeof p.abundance === 'number'? p.abundance : 0)));
      const gamma = 0.8;

      prods.forEach(p => {
        const yy = yForMW(p.mw);
        const norm = Math.pow(Math.max(0, p.abundance || 0) / maxAb, gamma);
        const alpha = Math.min(1.0, 0.08 + 0.92 * norm);
        const bandW = Math.floor(laneWidth * 0.7);
        const bandH = 10; // fixed band thickness for clarity in figures

        // Soft band with slight glow
        const x1 = cx - bandW/2;
        ctx.globalAlpha = alpha;
        const grad = ctx.createLinearGradient(0, yy-8, 0, yy+8);
        grad.addColorStop(0, 'rgba(255,255,255,0.12)');
        grad.addColorStop(0.5, 'rgba(255,255,255,1.0)');
        grad.addColorStop(1, 'rgba(255,255,255,0.12)');
        ctx.fillStyle = grad;
        ctx.shadowBlur = 2 + (norm * 8);
        ctx.shadowColor = 'white';
        ctx.fillRect(x1, yy - bandH/2, bandW, bandH);
        ctx.shadowBlur = 0;
        ctx.globalAlpha = 1;
      });

      // Labels under lane
      const top = lane.experiment || 'Experiment';
      let bottom = lane.name || '';
      ctx.save();
      ctx.translate(cx - 24, y + h - 16);
      ctx.rotate(-Math.PI/4);
      ctx.fillStyle = 'rgba(0,0,0,0.6)';
      ctx.fillRect(-6, -22, 120, 26);
      ctx.fillStyle = '#e5e7eb';
      ctx.font = '11px sans-serif';
      ctx.fillText(top.length > 16 ? top.slice(0,15)+'…' : top, 0, -6);
      ctx.fillText(bottom.length > 16 ? bottom.slice(0,15)+'…' : bottom, 0, 8);
      ctx.restore();
    });

    ctx.restore();
  }

  // Simple construct map (colored blocks)
  function drawConstructMap(ctx, x, y, w, h, assembled){
    if (!assembled || !assembled.sequence) return;
    const regions = assembled.regions || [];
    ctx.save();
    ctx.strokeStyle = '#cbd5e1';
    ctx.strokeRect(x, y, w, h);
    const total = assembled.sequence.length;
    const barY = y + px(h/3);
    const barH = px(h/3);
    const colorFor = (type) => {
      if (!type) return '#64748b';
      if (String(type).startsWith('RLUC')) return '#2563EB';
      switch(type){
        case 'FLUC': return '#DC2626';
        case 'LINKER': return '#F59E0B';
        case '5UTR': return '#94A3B8';
        case '3UTR': return '#CBD5E1';
        case 'CUSTOM': return '#6B7280';
        default: return '#64748b';
      }
    };
    regions.forEach(r => {
      const rs = (r.assembled ? r.start : (r.start - 1));
      const re = r.end;
      const x1 = x + Math.floor((rs/total) * w);
      const x2 = x + Math.floor((re/total) * w);
      const cw = Math.max(1, x2 - x1);
      const col = r.color || colorFor(r.type);
      ctx.fillStyle = col;
      ctx.globalAlpha = 0.2;
      ctx.fillRect(x1, barY, cw, barH);
      ctx.globalAlpha = 1;
      ctx.strokeStyle = col;
      ctx.lineWidth = 2;
      ctx.strokeRect(x1, barY, cw, barH);
      // label
      ctx.fillStyle = col;
      ctx.font = 'bold 11px sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText(String(r.type||'REGION'), x1 + cw/2, barY + barH + 14);
    });
    ctx.restore();
  }

  // Metrics list for figure
  function computeMetrics(lane, antibody){
    const products = Array.isArray(lane.products) ? lane.products : [];
    // Aggregate by reporters
    const rluc = products.filter(p => p.reporters && p.reporters.includes('RLUC')).reduce((a,b)=>a+(b.abundance||0),0);
    const fluc = products.filter(p => p.reporters && p.reporters.includes('FLUC')).reduce((a,b)=>a+(b.abundance||0),0);
    const bands = products.length;
    const translons = Array.isArray(lane.translons) ? lane.translons : [];
    const topT = [...translons].sort((a,b)=> (b.predictedAbundance||0)-(a.predictedAbundance||0)).slice(0,3);
    return { rluc, fluc, bands, topT };
  }

  function drawMetricsPanel(ctx, x, y, w, h, lane, antibody){
    ctx.save();
    ctx.strokeStyle = '#cbd5e1';
    ctx.strokeRect(x, y, w, h);
    const { rluc, fluc, bands, topT } = computeMetrics(lane, antibody);
    const lines = [
      `Construct: ${lane.experiment||'Experiment'} – ${lane.name||'Current'}`,
      `RLUC total: ${(rluc*100).toFixed(1)}%`,
      `FLUC total: ${(fluc*100).toFixed(1)}%`,
      `Predicted bands: ${bands}`,
      `Top translons:`
    ];
    ctx.fillStyle = '#0f172a';
    ctx.font = '12px sans-serif';
    let yy = y + 22;
    lines.forEach(s => { ctx.fillText(s, x + 12, yy); yy += 18; });
    topT.forEach(t => {
      const mw = ((t.endNt - t.startNt)/3*0.110).toFixed(1);
      const p = ((t.predictedAbundance||t.probability||0)*100).toFixed(1);
      const label = `${t.name||''}  start ${t.startNt}  frame ${t.frame}  ${mw} kDa  ${p}%`;
      ctx.fillText(label, x + 20, yy); yy += 16;
    });
    ctx.restore();
  }

  // Draw a panel label like "A", "B" in a small box
  function drawPanelLabel(ctx, x, y, letter){
    ctx.save();
    ctx.fillStyle = '#111827';
    ctx.fillRect(x, y, 26, 26);
    ctx.fillStyle = '#fff';
    ctx.font = 'bold 16px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText(letter, x+13, y+13);
    ctx.restore();
  }

  // Main API: compose + download PNG
  function exportFigurePNG(opts){
    const lanes = Array.isArray(opts.lanes) ? opts.lanes : [];
    if (lanes.length === 0) return;
    const focusIndex = Math.max(0, Math.min(lanes.length-1, (opts.focusIndex||lanes.length-1)));
    const focus = lanes[focusIndex];
    const antibody = opts.antibody || 'BOTH';

    // Figure size (hi-res)
    const W = 2200, H = 1400; // px
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');

    // White background
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(0, 0, W, H);
    ctx.fillStyle = '#0f172a';
    ctx.font = '18px sans-serif';
    ctx.fillText('Multi-panel Translation Figure', 24, 30);

    // Layout rects
    const pad = 24;
    const col1 = pad, col1W = Math.floor(W*0.62) - pad*1.5; // RDG big
    const col2 = Math.floor(W*0.62), col2W = W - col2 - pad; // gel
    const row1 = 50, row1H = Math.floor(H*0.56) - row1; // top row
    const row2 = Math.floor(H*0.56), row2H = H - row2 - pad; // bottom row

    // Panel A: RDG overview (offscreen render for the focused lane)
    if (root.RDGViz && focus && focus.assembled && focus.assembled.sequence) {
      const off = document.createElement('canvas');
      off.width = col1W; off.height = row1H;
      const offCtx = off.getContext('2d');
      try {
        const seq = focus.assembled.sequence;
        const feats = (root.RDGEngine && root.RDGEngine.identifyFeatures) ? root.RDGEngine.identifyFeatures(seq) : { predicted: { startCodons: [] } };
        const starts = feats.predicted.startCodons || (root.RDGEngine && root.RDGEngine.findStartCodons ? root.RDGEngine.findStartCodons(seq) : []);
        const trans = focus.translons || [];
        root.RDGViz.drawTreeLayout(offCtx, off, trans, seq, starts, ['#10b981','#3b82f6','#f59e0b'], focus.assembled.regions, (root.wbReadthroughStops||[]), []);
        ctx.drawImage(off, col1, row1);
      } catch (e) {
        offCtx.fillStyle = '#ef4444'; offCtx.fillText('RDG render failed', 10, 20); ctx.drawImage(off, col1, row1);
      }
      drawPanelLabel(ctx, col1, row1 - 22, 'A');
      ctx.font = '14px sans-serif'; ctx.fillStyle = '#334155'; ctx.fillText('RDG overview', col1 + 34, row1 - 6);
    }

    // Panel B: Gel lanes (all constructs)
    drawGelPanel(ctx, col2, row1, col2W, row1H, lanes);
    drawPanelLabel(ctx, col2, row1 - 22, 'B');
    ctx.font = '14px sans-serif'; ctx.fillStyle = '#334155'; ctx.fillText('Predicted Western blot', col2 + 34, row1 - 6);

    // Panel C: Construct map (focused)
    drawConstructMap(ctx, col1, row2, Math.floor(col1W*0.58), Math.floor(row2H*0.55), focus.assembled);
    drawPanelLabel(ctx, col1, row2 - 22, 'C');
    ctx.font = '14px sans-serif'; ctx.fillStyle = '#334155'; ctx.fillText('Construct map', col1 + 34, row2 - 6);

    // Panel D: Metrics (focused)
    const metricsX = col1 + Math.floor(col1W*0.60);
    drawMetricsPanel(ctx, metricsX, row2, col1 + col1W - metricsX, Math.floor(row2H*0.55), focus, antibody);
    drawPanelLabel(ctx, metricsX, row2 - 22, 'D');
    ctx.font = '14px sans-serif'; ctx.fillStyle = '#334155'; ctx.fillText('Summary metrics', metricsX + 34, row2 - 6);

    // Caption area (bottom)
    const capY = row2 + Math.floor(row2H*0.62);
    ctx.fillStyle = '#0f172a'; ctx.font = '12px sans-serif';
    const cap = opts.caption || 'RDG shows translation paths; gel summarizes predicted detectable products; metrics list dominant translons and reporter totals.';
    ctx.fillText(cap, col1, capY);

    // Download
    const a = document.createElement('a');
    a.download = opts.filename || 'translation_figure.png';
    a.href = canvas.toDataURL('image/png');
    a.click();
  }

  const api = { exportFigurePNG };
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = api;
  } else {
    root.RDGFigure = api;
  }
})(typeof window !== 'undefined' ? window : global);

