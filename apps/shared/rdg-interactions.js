/**
 * RDG Interactions - shared click handling for RDG canvases (non-simulation)
 *
 * Expects RDG viz to expose on window:
 *  - RDG_LAST_EDGES: array of edges with type ('translation', 'readthrough', ...)
 *  - RDG_NT_TO_PIXEL: function(nt)->x
 *  - RDG_FRAME_Y: [y0,y1,y2]
 *  - RDG_Y_OFFSET: number
 */
(function (root) {
  function hitTestTranslonAt(x, y) {
    const edges = root.RDG_LAST_EDGES || [];
    const ntToPx = root.RDG_NT_TO_PIXEL || ((nt) => nt);
    const yOff = root.RDG_Y_OFFSET || 0;
    const tol = 40; // generous tolerance to ensure clicks register
    for (const e of edges) {
      if (e.type !== 'translation') continue;
      const x1 = ntToPx(e.x1);
      const x2 = ntToPx(e.x2);
      const yLine = (e.y1 || 0) + yOff;
      if (x >= x1 && x <= x2 && Math.abs(y - yLine) <= tol) {
        return e.translon || null;
      }
    }
    return null;
  }

  function hitTestStopAt(x, y, sequence) {
    if (!sequence) return null;
    const rails = root.RDG_FRAME_Y || [80, 110, 140];
    const ntToPx = root.RDG_NT_TO_PIXEL || ((nt) => nt);
    // Determine clicked frame by nearest rail within 15px
    let frame = -1, minDy = 1e9;
    rails.forEach((yy, idx) => {
      const d = Math.abs(y - yy);
      if (d < minDy && d < 15) { minDy = d; frame = idx; }
    });
    if (frame === -1) return null;
    // Map x to nearest nt
    // Invert ntToPx using a coarse search (sequence length can be large; do local vicinity)
    // We'll scan stops in this frame and choose the one with minimal |px - x|
    const stops = [];
    for (let pos = frame; pos < sequence.length - 2; pos += 3) {
      const codon = sequence.substring(pos, pos + 3);
      if (codon === 'UAA' || codon === 'UAG' || codon === 'UGA') {
        stops.push({ pos, codon, frame });
      }
    }
    let best = null, bestDx = 1e9;
    for (const s of stops) {
      const px = ntToPx(s.pos);
      const dx = Math.abs(px - x);
      if (dx < bestDx && dx < 12) { bestDx = dx; best = s; }
    }
    return best;
  }

  // Hit-test near start codon tick markers
  function hitTestStartAt(x, y) {
    const starts = root.RDG_LAST_STARTS || [];
    const rails = root.RDG_FRAME_Y || [80, 110, 140];
    const ntToPx = root.RDG_NT_TO_PIXEL || ((nt) => nt);
    // Find closest start within 10px in X and 12px to its frame Y
    let best = null, bestDx = 1e9;
    starts.forEach(s => {
      const sx = ntToPx(s.pos);
      const sy = rails[s.frame] || 0;
      const dx = Math.abs(x - sx);
      const dy = Math.abs(y - sy);
      if (dx <= 10 && dy <= 12 && dx < bestDx) { bestDx = dx; best = s; }
    });
    return best; // {pos, frame, codon, ...}
  }

  function wireInteractions(canvas, opts) {
    if (!canvas) return;
    const getMode = opts.getMode || (() => 'translon');
    const getSequence = opts.getSequence || (() => '');
    const onAdd = opts.onAddTranslon || (() => {});
    const onEdit = opts.onEditTranslon || (() => {});
    const onToggleRT = opts.onToggleReadthrough || (() => {});

    canvas.addEventListener('click', (e) => {
      const rect = canvas.getBoundingClientRect();
      const scaleX = canvas.width / rect.width;
      const scaleY = canvas.height / rect.height;
      const x = (e.clientX - rect.left) * scaleX;
      const y = (e.clientY - rect.top) * scaleY;
      const mode = getMode();

      if (mode === 'edit') {
        const t = hitTestTranslonAt(x, y);
        if (t) onEdit(t);
        return;
      }

      if (mode === 'readthrough') {
        const seq = getSequence();
        const s = hitTestStopAt(x, y, seq);
        if (s) onToggleRT(s);
        return;
      }

      // default: add translon near a start on clicked frame
      if (typeof opts.onAddTranslonAt !== 'function') return;
      opts.onAddTranslonAt(x, y);
    });
  }

  const api = { hitTestTranslonAt, hitTestStopAt, hitTestStartAt, wireInteractions };
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = api;
  } else {
    root.RDGInteractions = api;
  }
})(typeof window !== 'undefined' ? window : global);
