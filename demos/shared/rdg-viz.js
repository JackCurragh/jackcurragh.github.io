/**
 * RDG Visualization Module
 *
 * Provides tree and flux visualization rendering for RDG analysis
 * Shared between RDG demo and Western Blot demo
 */

// Calculate tree layout from translons with support for readthrough and frameshifts
function calculateTreeLayout(translons, sequence, readthroughStops = [], frameshiftSites = []) {
    console.log(`[calculateTreeLayout] Called with ${translons.length} translons, ${readthroughStops.length} readthrough stops, ${frameshiftSites.length} frameshift sites`);
    console.log(`[calculateTreeLayout] Readthrough stops:`, readthroughStops);

    const nodes = [];
    const edges = [];
    const branchSpacing = 80;

    let currentY = 200;
    const rootNode = {
        id: 'root',
        x: 0,
        y: currentY,
        type: 'root'
    };
    nodes.push(rootNode);

    let previousNode = rootNode;

    if (translons.length === 0) {
        const endNode = {
            id: 'end',
            x: sequence.length,
            y: currentY,
            type: 'end'
        };
        nodes.push(endNode);

        edges.push({
            from: 'root',
            to: 'end',
            x1: 0,
            y1: currentY,
            x2: sequence.length,
            y2: currentY,
            type: 'noncoding'
        });

        return { nodes, edges };
    }

    const sortedTranslons = [...translons].sort((a, b) => a.startNt - b.startNt);

    // Helper function to calculate the Y position of the non-coding line at a given X position
    function getScanningYAtPosition(xPos) {
        // Start at initial Y
        let yPos = 200;

        // For each translon whose start is BEFORE xPos, add branchSpacing/2
        for (const t of sortedTranslons) {
            if (t.startNt <= xPos) {
                yPos += branchSpacing / 2;
            } else {
                break; // translons are sorted, so we can stop
            }
        }

        return yPos;
    }

    sortedTranslons.forEach((translon, index) => {
        const branchX = translon.startNt;
        const branchY = currentY;

        const decisionNode = {
            id: `decision_${index}`,
            x: branchX,
            y: branchY,
            type: 'decision'
        };
        nodes.push(decisionNode);

        edges.push({
            from: previousNode.id,
            to: decisionNode.id,
            x1: previousNode.x,
            y1: previousNode.y,
            x2: branchX,
            y2: branchY,
            type: 'noncoding'
        });

        // Check if this translon has a readthrough stop within it
        console.log(`Checking translon ${translon.name}: startNt=${translon.startNt}, endNt=${translon.endNt}, frame=${translon.frame}`);
        console.log(`  Available readthrough stops:`, readthroughStops);

        const readthroughInTranslon = readthroughStops.find(rt => {
            const isInRange = rt.pos >= translon.startNt && rt.pos < translon.endNt;
            const isInFrame = rt.frame === translon.frame;
            console.log(`    Checking RT at pos=${rt.pos}, frame=${rt.frame}: inRange=${isInRange}, inFrame=${isInFrame}`);
            return isInRange && isInFrame;
        });

        console.log(`  Result: readthrough=${readthroughInTranslon ? readthroughInTranslon.pos : 'NONE'}`);

        if (readthroughInTranslon) {
            console.log(`  ✓✓✓ Creating readthrough translon - marking stop at position ${readthroughInTranslon.pos}`);

            // Find the NEXT in-frame stop AFTER the readthrough stop
            let nextStopPos = sequence.length;
            const frame = translon.frame;
            for (let pos = readthroughInTranslon.pos + 3; pos < sequence.length - 2; pos += 3) {
                const codon = sequence.substring(pos, pos + 3);
                if (['UAA', 'UAG', 'UGA'].includes(codon)) {
                    // Check if THIS stop is also marked for readthrough
                    const isAlsoReadthrough = readthroughStops.some(rt => rt.pos === pos);
                    if (!isAlsoReadthrough) {
                        nextStopPos = pos + 3; // Include the stop codon
                        break;
                    }
                }
            }

            console.log(`  Extended translon to next stop at ${nextStopPos} (was ${translon.endNt})`);

            // Readthrough case: draw extended translon and mark where the first stop is
            const translationY = branchY - branchSpacing / 2;
            const translationEndNode = {
                id: `${translon.name}_end`,
                x: nextStopPos,
                y: translationY,
                type: 'endpoint',
                translon: translon
            };
            nodes.push(translationEndNode);

            edges.push({
                from: decisionNode.id,
                to: translationEndNode.id,
                x1: branchX,
                y1: branchY,
                x2: branchX,
                y2: translationY,
                type: 'vertical_branch'
            });

            // Draw the full extended translation bar with readthrough marker
            edges.push({
                from: decisionNode.id,
                to: translationEndNode.id,
                x1: branchX,
                y1: translationY,
                x2: nextStopPos,
                y2: translationY,
                type: 'translation',
                translon: translon,
                readthroughStop: readthroughInTranslon.pos, // Mark where the first stop is
                readthroughProb: readthroughInTranslon.probability
            });

            // Scanning continues (goes down)
            const scanningY = branchY + branchSpacing / 2;
            const scanningNode = {
                id: `scanning_${index}`,
                x: branchX,
                y: scanningY,
                type: 'scanning'
            };
            nodes.push(scanningNode);

            edges.push({
                from: decisionNode.id,
                to: scanningNode.id,
                x1: branchX,
                y1: branchY,
                x2: branchX,
                y2: scanningY,
                type: 'vertical_branch'
            });

            // Reinitiation path 1: From first stop straight down to non-coding line at that X position
            const firstStopPos = readthroughInTranslon.pos + 3;
            const firstStopScanningY = getScanningYAtPosition(firstStopPos);
            edges.push({
                from: `${translon.name}_first_stop`,
                to: `${translon.name}_first_stop_scan`,
                x1: firstStopPos,
                y1: translationY,
                x2: firstStopPos,
                y2: firstStopScanningY,
                type: 'reinitiation',
                translon: translon
            });

            // Reinitiation path 2: From extended end straight down to non-coding line at that X position
            const extendedStopScanningY = getScanningYAtPosition(nextStopPos);
            edges.push({
                from: translationEndNode.id,
                to: `${translon.name}_extended_stop_scan`,
                x1: nextStopPos,
                y1: translationY,
                x2: nextStopPos,
                y2: extendedStopScanningY,
                type: 'reinitiation',
                translon: translon
            });

            previousNode = scanningNode;
            currentY = scanningY;

        } else {
            // Normal case: no readthrough
            const translationY = branchY - branchSpacing / 2;
            const translationEndNode = {
                id: `${translon.name}_end`,
                x: translon.endNt,
                y: translationY,
                type: 'endpoint',
                translon: translon
            };
            nodes.push(translationEndNode);

            edges.push({
                from: decisionNode.id,
                to: translationEndNode.id,
                x1: branchX,
                y1: branchY,
                x2: branchX,
                y2: translationY,
                type: 'vertical_branch'
            });

            edges.push({
                from: decisionNode.id,
                to: translationEndNode.id,
                x1: branchX,
                y1: translationY,
                x2: translon.endNt,
                y2: translationY,
                type: 'translation',
                translon: translon
            });

            // Scanning continues (goes down)
            const scanningY = branchY + branchSpacing / 2;
            const scanningNode = {
                id: `scanning_${index}`,
                x: branchX,
                y: scanningY,
                type: 'scanning'
            };
            nodes.push(scanningNode);

            edges.push({
                from: decisionNode.id,
                to: scanningNode.id,
                x1: branchX,
                y1: branchY,
                x2: branchX,
                y2: scanningY,
                type: 'vertical_branch'
            });

            // Reinitiation path - straight down to non-coding line at stop position
            const stopScanningY = getScanningYAtPosition(translon.endNt);
            edges.push({
                from: translationEndNode.id,
                to: scanningNode.id,
                x1: translon.endNt,
                y1: translationY,
                x2: translon.endNt,
                y2: stopScanningY,
                type: 'reinitiation',
                translon: translon
            });

            previousNode = scanningNode;
            currentY = scanningY;
        }
    });

    const endNode = {
        id: 'end',
        x: sequence.length,
        y: currentY,
        type: 'end'
    };
    nodes.push(endNode);

    edges.push({
        from: previousNode.id,
        to: endNode.id,
        x1: previousNode.x,
        y1: previousNode.y,
        x2: sequence.length,
        y2: currentY,
        type: 'noncoding'
    });

    return { nodes, edges };
}

// Draw tree layout on canvas
function drawTreeLayout(ctx, canvas, translons, sequence, startCodons, FRAME_COLORS, constructRegions = null, readthroughStops = [], frameshiftSites = []) {
    console.log(`[drawTreeLayout] VERSION 2.0 - Called with ${readthroughStops.length} readthrough stops`);
    console.log(`[drawTreeLayout] Readthrough stops:`, readthroughStops);

    // Config (supports accessibility/color modes)
    const getConfig = () => (typeof window !== 'undefined' && window.RDGVizConfig) ? window.RDGVizConfig : { mode: 'accessibility' };
    const cfg = getConfig();

    // Palette (construct blocks)
    const colorForType = (type) => {
        if (!type) return '#64748b';
        if (typeof type === 'string' && type.startsWith('RLUC')) return '#2563EB'; // blue-600
        switch (type) {
            case 'FLUC': return '#DC2626'; // red-600
            case 'LINKER': return '#F59E0B'; // amber-500
            case '5UTR': return '#94A3B8'; // slate-300 (context)
            case '3UTR': return '#CBD5E1'; // slate-200 (context)
            case 'CUSTOM': return '#6B7280'; // grey-500
            default: return '#64748b';
        }
    };
    const baseTranslonColor = cfg.mode === 'frames' ? null : '#475569'; // slate-600
    const railFill = cfg.mode === 'frames' ? null : '#E5E7EB'; // gray-200
    const railStroke = cfg.mode === 'frames' ? null : '#CBD5E1'; // gray-300
    const showOrfBarsOnFrames = !!cfg.showOrfBarsOnFrames;

    const margin = 50;
    const seqEndX = canvas.width - margin;
    const ntToPixel = (nt) => margin + (nt / sequence.length) * (seqEndX - margin);

    console.log(`[drawTreeLayout] About to call calculateTreeLayout...`);
    // Call the internal calculateTreeLayout function directly (not from global scope)
    const result = (function() {
        return calculateTreeLayout(translons, sequence, readthroughStops, frameshiftSites);
    })();
    const { nodes, edges } = result;
    console.log(`[drawTreeLayout] calculateTreeLayout returned ${nodes.length} nodes and ${edges.length} edges`);

    // Ensure canvas is tall enough for all nodes/edges
    let maxY = 0;
    nodes.forEach(n => { if (n.y > maxY) maxY = n.y; });
    edges.forEach(e => { maxY = Math.max(maxY, e.y1 || 0, e.y2 || 0); });
    const neededH = Math.max(400, Math.ceil(maxY + 80));
    if (canvas.height < neededH) {
        canvas.height = neededH;
    }

    // Frame rail Y positions (padding under construct bar)
    const railsY = [90, 120, 150];

    // Highlight sets (global, optional)
    const hiKeys = (typeof window !== 'undefined' && window.RDG_HIGHLIGHT_KEYS) ? window.RDG_HIGHLIGHT_KEYS : null;
    const hiReporters = (typeof window !== 'undefined' && window.RDG_HIGHLIGHT_REPORTERS) ? window.RDG_HIGHLIGHT_REPORTERS : null;

    // Draw construct regions if provided
    if (constructRegions && constructRegions.length > 0) {
        constructRegions.forEach(region => {
            const rStartNt = region.assembled ? region.start : (region.start - 1);
            const rEndNt = region.end;
            const x1 = ntToPixel(rStartNt);
            const x2 = ntToPixel(rEndNt);
            const regionColor = region.color || colorForType(region.type);

            // Draw on top bar
            ctx.fillStyle = regionColor;
            ctx.globalAlpha = 0.15;
            ctx.fillRect(x1, 20, x2 - x1, 40);
            ctx.globalAlpha = 1;

            ctx.strokeStyle = regionColor;
            ctx.lineWidth = 2;
            ctx.strokeRect(x1, 20, x2 - x1, 40);

            // Highlight entire reporter block if requested
            if (hiReporters && hiReporters.has(region.type)) {
                ctx.save();
                ctx.strokeStyle = '#22c55e';
                ctx.lineWidth = 4;
                ctx.globalAlpha = 0.9;
                ctx.strokeRect(x1 - 2, 18, (x2 - x1) + 4, 44);
                ctx.restore();
            }

            ctx.fillStyle = regionColor;
            ctx.font = 'bold 10px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(region.type, (x1 + x2) / 2, 50);

            // Optional: ORF bar on frame track (off by default in accessibility mode)
            if (showOrfBarsOnFrames && (region.type === 'RLUC' || region.type === 'RLUC_WEAK' || region.type === 'RLUC_NO_STOP' || region.type === 'FLUC')) {
                const frameIdx = (rStartNt) % 3;
                const y = railsY[frameIdx];
                ctx.fillStyle = regionColor;
                ctx.globalAlpha = 0.35;
                ctx.fillRect(x1, y - 10, x2 - x1, 20);
                ctx.globalAlpha = 1;
                ctx.strokeStyle = regionColor;
                ctx.lineWidth = 3;
                ctx.strokeRect(x1, y - 10, x2 - x1, 20);
            }
        });
    }

    // Draw frame tracks (neutral in accessibility mode)
    const frameY = railsY;
    const frameBarHeight = 20;
    const frameLabels = ['Frame 0', 'Frame 1', 'Frame 2'];

    ctx.font = 'bold 11px sans-serif';
    ctx.textAlign = 'right';

    frameY.forEach((y, frame) => {
        const railFillColor = cfg.mode === 'frames' ? FRAME_COLORS[frame] : railFill;
        const railStrokeColor = cfg.mode === 'frames' ? FRAME_COLORS[frame] : railStroke;
        const railAlpha = cfg.mode === 'frames' ? 0.2 : 1.0;

        ctx.fillStyle = '#334155';
        ctx.fillText(frameLabels[frame], margin - 10, y + 6);

        if (railFillColor) {
            ctx.fillStyle = railFillColor;
            ctx.globalAlpha = railAlpha;
            ctx.fillRect(margin, y - frameBarHeight/2, seqEndX - margin, frameBarHeight);
            ctx.globalAlpha = 1;
        }
        if (railStrokeColor) {
            ctx.strokeStyle = railStrokeColor;
            ctx.lineWidth = 1;
            ctx.strokeRect(margin, y - frameBarHeight/2, seqEndX - margin, frameBarHeight);
        }
    });

    // Expose last layout and helpers for interactivity (only for primary RDG canvas)
    if (typeof window !== 'undefined') {
        const isPrimary = canvas && (canvas.id === 'rdg-canvas' || window.__RDG_PRIMARY_CANVAS === canvas);
        if (isPrimary) {
            window.RDG_LAST_EDGES = edges;
            window.RDG_NT_TO_PIXEL = ntToPixel;
            window.RDG_FRAME_Y = frameY;
            window.RDG_Y_OFFSET = (typeof rdgYOffset !== 'undefined') ? rdgYOffset : 0;
            window.RDG_LAST_STARTS = startCodons || [];
        }
    }

    // Draw start codons
    for (const start of startCodons) {
        const x = ntToPixel(start.pos);
        const y = frameY[start.frame];

        if (start.isAUG) {
            ctx.strokeStyle = '#059669';
            ctx.lineWidth = 3;
        } else {
            ctx.strokeStyle = '#34d399';
            ctx.lineWidth = 2;
        }
        ctx.beginPath();
        ctx.moveTo(x, y - frameBarHeight/2);
        ctx.lineTo(x, y + frameBarHeight/2);
        ctx.stroke();
    }

    // Draw ALL stop codons in all frames (not just those from startCodons)
    for (let frame = 0; frame < 3; frame++) {
        for (let pos = frame; pos < sequence.length - 2; pos += 3) {
            const codon = sequence.substring(pos, pos + 3);
            if (['UAA', 'UAG', 'UGA'].includes(codon)) {
                const stopX = ntToPixel(pos);
                const y = frameY[frame];
                ctx.strokeStyle = '#dc2626';
                ctx.lineWidth = 2;
                ctx.beginPath();
                ctx.moveTo(stopX, y - frameBarHeight/2);
                ctx.lineTo(stopX, y + frameBarHeight/2);
                ctx.stroke();
            }
        }
    }

    // Draw readthrough stop codons (with special marker)
    readthroughStops.forEach(rt => {
        const x = ntToPixel(rt.pos);
        const y = frameY[rt.frame];

        // Draw as hollow circle to indicate potential readthrough
        ctx.strokeStyle = '#f59e0b'; // Orange for readthrough
        ctx.lineWidth = 3;
        ctx.beginPath();
        ctx.arc(x, y, 8, 0, Math.PI * 2);
        ctx.stroke();

        // Draw vertical line through it
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(x, y - frameBarHeight/2);
        ctx.lineTo(x, y + frameBarHeight/2);
        ctx.stroke();

        // Add "RT" label
        ctx.fillStyle = '#f59e0b';
        ctx.font = 'bold 10px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('RT', x, y - 15);
    });

    // Draw frameshift sites (with arrow indicator)
    frameshiftSites.forEach(fs => {
        const x = ntToPixel(fs.pos);
        const fromY = frameY[fs.fromFrame];
        const toY = frameY[fs.toFrame];

        // Draw arrow from one frame to another
        ctx.strokeStyle = '#8b5cf6'; // Purple for frameshift
        ctx.lineWidth = 3;
        ctx.setLineDash([5, 3]);
        ctx.beginPath();
        ctx.moveTo(x, fromY);
        ctx.lineTo(x, toY);
        ctx.stroke();
        ctx.setLineDash([]);

        // Arrow head
        ctx.fillStyle = '#8b5cf6';
        ctx.beginPath();
        if (fs.toFrame > fs.fromFrame) {
            // Downward arrow
            ctx.moveTo(x, toY);
            ctx.lineTo(x - 4, toY - 6);
            ctx.lineTo(x + 4, toY - 6);
        } else {
            // Upward arrow
            ctx.moveTo(x, toY);
            ctx.lineTo(x - 4, toY + 6);
            ctx.lineTo(x + 4, toY + 6);
        }
        ctx.closePath();
        ctx.fill();

        // Add "FS" label
        ctx.fillStyle = '#8b5cf6';
        ctx.font = 'bold 10px sans-serif';
        ctx.textAlign = 'center';
        const labelY = Math.min(fromY, toY) - 15;
        ctx.fillText(`FS${fs.shift > 0 ? '+' : ''}${fs.shift}`, x, labelY);
    });

    // Draw 5' and 3' labels
    ctx.fillStyle = '#64748b';
    ctx.font = 'bold 14px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText("5'", margin - 30, frameY[1] + 4);
    ctx.fillText("3'", seqEndX + 20, frameY[1] + 4);

    // Draw RDG edges
    edges.forEach(edge => {
        const x1 = ntToPixel(edge.x1);
        const x2 = ntToPixel(edge.x2);

        if (edge.type === 'noncoding') {
            const y1p = (edge.y1 || 0) + (window.RDG_Y_OFFSET || 0);
            const y2p = (edge.y2 || 0) + (window.RDG_Y_OFFSET || 0);
            ctx.strokeStyle = '#999';
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(x1, y1p);
            ctx.lineTo(x2, y2p);
            ctx.stroke();
        } else if (edge.type === 'vertical_branch') {
            ctx.strokeStyle = '#000';
            ctx.lineWidth = 3;
            ctx.beginPath();
            ctx.moveTo(x1, edge.y1);
            ctx.lineTo(x2, edge.y2);
            ctx.stroke();
        } else if (edge.type === 'translation') {
            // Base: translon bar (apply RDG offset)
            const yOff = (typeof window !== 'undefined' && typeof window.RDG_Y_OFFSET === 'number') ? window.RDG_Y_OFFSET : 0;
            ctx.strokeStyle = (cfg.mode === 'frames' ? FRAME_COLORS[edge.translon.frame] : baseTranslonColor);
            ctx.lineWidth = 12;
            ctx.lineCap = 'butt';
            ctx.beginPath();
            ctx.moveTo(x1, (edge.y1||0) + yOff);
            ctx.lineTo(x2, (edge.y2||0) + yOff);
            ctx.stroke();

            // (moved) highlight overlay is drawn after overlays/labels so it stays on top

            // NOTE: Do not overlay in-situ sequence (UTRs/CUSTOM). Overlay for reporter/linker only (handled below).

            // Outlines: indicate construct blocks (reporters/linker) without obscuring frame color
            if (constructRegions && constructRegions.length) {
                const overlayTypes = new Set(['RLUC','RLUC_WEAK','RLUC_NO_STOP','FLUC','LINKER']);
                constructRegions.forEach(region => {
                    if (!overlayTypes.has(region.type)) return;
                    const rStartNt = region.assembled ? region.start : (region.start - 1);
                    const rEndNt = region.end;
                    const segStartNt = Math.max(edge.x1, rStartNt);
                    const segEndNt = Math.min(edge.x2, rEndNt);
                    if (segEndNt > segStartNt) {
                        const sx = ntToPixel(segStartNt);
                        const ex = ntToPixel(segEndNt);
                        const color = region.color || colorForType(region.type);
                        const barHalf = 8;
                        const yTop = edge.y1 - barHalf;
                        const height = barHalf * 2;
                        // Strongly mute base within span
                        ctx.save();
                        ctx.globalAlpha = 0.7;
                        ctx.fillStyle = '#ffffff';
                        ctx.fillRect(sx, yTop + 2, ex - sx, height - 4);
                        ctx.restore();
                        // Thicker outline
                        ctx.strokeStyle = color;
                        ctx.lineWidth = 5;
                        ctx.setLineDash([]);
                        ctx.strokeRect(sx, yTop, ex - sx, height);
                    }
                });
            }

            // Determine which reporters this translon produces
            let reporters = [];
            if (constructRegions) {
                constructRegions.forEach(region => {
                    if ((region.type === 'RLUC' || region.type === 'RLUC_WEAK' || region.type === 'RLUC_NO_STOP' || region.type === 'FLUC') &&
                        edge.translon.startNt <= ((region.assembled ? region.start : (region.start - 1))) &&
                        edge.translon.endNt >= region.end) {
                        reporters.push(region.type);
                    }
                });
            }

            // Label with reporter info
            const displayName = edge.translon.displayName || edge.translon.name;
            ctx.fillStyle = FRAME_COLORS[edge.translon.frame];
            ctx.font = 'bold 11px sans-serif';
            ctx.textAlign = 'left';
            const labelX = x1 + 5;
            const labelY = ((edge.y1||0) + ((typeof window !== 'undefined' && typeof window.RDG_Y_OFFSET === 'number') ? window.RDG_Y_OFFSET : 0)) - 8;

            let labelText = `${displayName}`;
            if (reporters.length > 0) {
                labelText += ` → ${reporters.join('+')}`;
            }
            ctx.fillText(labelText, labelX, labelY);

            // Draw readthrough stop marker if this translon has readthrough
            if (edge.readthroughStop !== undefined) {
                const stopX = ntToPixel(edge.readthroughStop);
                // Draw vertical line at the readthrough stop position
                ctx.strokeStyle = '#f59e0b';
                ctx.lineWidth = 3;
                ctx.setLineDash([]);
                ctx.beginPath();
                const yOff2 = (typeof window !== 'undefined' && typeof window.RDG_Y_OFFSET === 'number') ? window.RDG_Y_OFFSET : 0;
                ctx.moveTo(stopX, ((edge.y1||0)+yOff2) - 6);
                ctx.lineTo(stopX, ((edge.y1||0)+yOff2) + 6);
                ctx.stroke();

                // Add "RT" label above
                ctx.fillStyle = '#f59e0b';
                ctx.font = 'bold 10px sans-serif';
                ctx.textAlign = 'center';
                ctx.fillText(`RT`, stopX, edge.y1 - 12);
            }

            // Draw connection lines to reporter regions
            if (reporters.length > 0 && constructRegions) {
                constructRegions.forEach(region => {
                    if (reporters.includes(region.type)) {
                        const regX1 = ntToPixel(region.start - 1);
                        const regX2 = ntToPixel(region.end);
                        const regFrame = (region.start - 1) % 3;
                        const regY = (typeof railsY !== 'undefined' ? railsY[regFrame] : [80,110,140][regFrame]);

                        // Dashed connection line from translation bar to reporter
                        ctx.strokeStyle = region.color;
                        ctx.lineWidth = 1;
                        ctx.setLineDash([3, 3]);
                        ctx.globalAlpha = 0.4;

                        // Connect from middle of translation bar to reporter
                        const midX = (Math.max(x1, regX1) + Math.min(x2, regX2)) / 2;
                        ctx.beginPath();
                        const y1p2 = (edge.y1 || 0) + (window.RDG_Y_OFFSET || 0);
                        ctx.moveTo(midX, y1p2 + 6);
                        ctx.lineTo(midX, regY);
                        ctx.stroke();

                        ctx.setLineDash([]);
                        ctx.globalAlpha = 1;
                    }
                });
            }

            // Highlight translon if selected (draw on top of overlays/labels)
            try {
                if (hiKeys) {
                    const k = `${edge.translon.startNt}|${edge.translon.endNt}|${edge.translon.frame}`;
                    if (hiKeys.has(k)) {
                        ctx.save();
                        ctx.strokeStyle = '#22c55e';
                        ctx.lineWidth = 14;
                        ctx.globalAlpha = 0.8;
                        ctx.beginPath();
                        ctx.moveTo(x1, (edge.y1||0) + yOff);
                        ctx.lineTo(x2, (edge.y2||0) + yOff);
                        ctx.stroke();
                        ctx.restore();
                    }
                }
            } catch {}
        } else if (edge.type === 'readthrough') {
            // Readthrough segment - draw with dotted line to show it's conditional
            ctx.strokeStyle = (cfg.mode === 'frames' ? FRAME_COLORS[edge.translon.frame] : baseTranslonColor);
            ctx.lineWidth = 12;
            ctx.lineCap = 'butt';
            ctx.setLineDash([5, 5]);
            ctx.globalAlpha = 0.7;
            ctx.beginPath();
            ctx.moveTo(x1, edge.y1);
            ctx.lineTo(x2, edge.y2);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.globalAlpha = 1;

            // Outlines for readthrough path (reporter/linker only)
            if (constructRegions && constructRegions.length) {
                const overlayTypes = new Set(['RLUC','RLUC_WEAK','RLUC_NO_STOP','FLUC','LINKER']);
                constructRegions.forEach(region => {
                    if (!overlayTypes.has(region.type)) return;
                    const rStartNt = region.assembled ? region.start : (region.start - 1);
                    const rEndNt = region.end;
                    const segStartNt = Math.max(edge.x1, rStartNt);
                    const segEndNt = Math.min(edge.x2, rEndNt);
                    if (segEndNt > segStartNt) {
                        const sx = ntToPixel(segStartNt);
                        const ex = ntToPixel(segEndNt);
                        const color = region.color || colorForType(region.type);
                        const barHalf = 8;
                        const yTop = edge.y1 - barHalf;
                        const height = barHalf * 2;
                        // Strongly mute base within span
                        ctx.save();
                        ctx.globalAlpha = 0.7;
                        ctx.fillStyle = '#ffffff';
                        ctx.fillRect(sx, yTop + 2, ex - sx, height - 4);
                        ctx.restore();
                        // Thicker dashed outline
                        ctx.strokeStyle = color;
                        ctx.lineWidth = 5;
                        ctx.setLineDash([5, 5]);
                        ctx.strokeRect(sx, yTop, ex - sx, height);
                        ctx.setLineDash([]);
                    }
                });
            }

            // Add readthrough probability label
            if (edge.probability !== undefined) {
                ctx.fillStyle = '#f59e0b';
                ctx.font = 'bold 10px sans-serif';
                ctx.textAlign = 'center';
                const midX = (x1 + x2) / 2;
                const y1p = (edge.y1 || 0) + (window.RDG_Y_OFFSET || 0);
                ctx.fillText(`RT ${(edge.probability * 100).toFixed(0)}%`, midX, y1p - 5);
            }
        } else if (edge.type === 'reinitiation') {
            const y1p = (edge.y1 || 0) + (window.RDG_Y_OFFSET || 0);
            const y2p = (edge.y2 || 0) + (window.RDG_Y_OFFSET || 0);
            console.log(`Drawing reinitiation: x1=${x1}, y1=${y1p}, x2=${x2}, y2=${y2p}, label=${edge.label}`);
            ctx.strokeStyle = '#8b5cf6';
            ctx.lineWidth = 2;
            ctx.setLineDash([8, 4]);
            ctx.globalAlpha = 0.6;
            ctx.beginPath();
            ctx.moveTo(x1, y1p);
            ctx.lineTo(x2, y2p);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.globalAlpha = 1;

            // Removed probability labels from reinitiation paths
        }
    });

    // Draw decision points
    nodes.forEach(node => {
        if (node.type === 'decision') {
            const x = ntToPixel(node.x);
            ctx.fillStyle = '#000';
            ctx.beginPath();
            const y = (node.y || 0) + (window.RDG_Y_OFFSET || 0);
            ctx.arc(x, y, 5, 0, Math.PI * 2);
            ctx.fill();
        } else if (node.type === 'readthrough_decision') {
            // Special marker for readthrough decision points
            const x = ntToPixel(node.x);
            ctx.strokeStyle = '#f59e0b';
            ctx.fillStyle = '#fef3c7';
            ctx.lineWidth = 3;
            ctx.beginPath();
            const y = (node.y || 0) + (window.RDG_Y_OFFSET || 0);
            ctx.arc(x, y, 7, 0, Math.PI * 2);
            ctx.fill();
            ctx.stroke();
        }
    });

    // Legend (optional)
    if (cfg.showLegend) {
        const legendX = canvas.width - 220;
        const legendY = 10;
        const drawChip = (x, y, w, h, color, label, dashed=false) => {
            if (dashed) { ctx.setLineDash([6,3]); }
            ctx.strokeStyle = color; ctx.lineWidth = 5; ctx.strokeRect(x, y, w, h);
            ctx.setLineDash([]);
            ctx.fillStyle = '#334155'; ctx.font = '10px sans-serif'; ctx.textAlign = 'left'; ctx.fillText(label, x + w + 6, y + h - 2);
        };
        ctx.fillStyle = 'rgba(255,255,255,0.9)'; ctx.strokeStyle = '#CBD5E1'; ctx.lineWidth = 1;
        ctx.fillRect(legendX, legendY, 210, 86); ctx.strokeRect(legendX, legendY, 210, 86);
        ctx.fillStyle = '#334155'; ctx.font = 'bold 11px sans-serif'; ctx.fillText('Legend', legendX + 8, legendY + 14);
        const lY = legendY + 26; const h = 10; const w = 34; let curX = legendX + 8;
        drawChip(curX, lY, w, h, colorForType('RLUC'), 'RLUC'); curX += 74;
        drawChip(curX, lY, w, h, colorForType('FLUC'), 'FLUC'); curX = legendX + 8; const lY2 = lY + 18;
        drawChip(curX, lY2, w, h, colorForType('LINKER'), 'LINKER'); curX += 92;
        drawChip(curX, lY2, w, h, colorForType('RLUC'), 'readthrough', true);
    }
}

// Suggest construct design based on sequence analysis
function suggestConstructDesign(sequence, startCodons, features = null) {
    const suggestions = [];

    // If features are available, use them to create targeted suggestions
    if (features && features.predicted) {
        const starts = features.predicted.startCodons || startCodons || [];

        // uORF-based construct: pick strongest uORF by initiation probability
        const uorfs = (features.predicted.uorfs || []).slice().sort((a, b) => (b.initiationProbability || 0) - (a.initiationProbability || 0));
        if (uorfs.length > 0) {
            const top = uorfs[0];
            const pad = 50;
            const lastUorf = top.end;
            const fiveEnd = Math.max(1, top.start - pad);
            const fiveStop = Math.min(sequence.length, lastUorf + pad);

            suggestions.push({
                type: 'uORF Reporter',
                reason: `Targets strongest uORF at ${top.start} (Kozak ${(top.kozakScore||0).toFixed(2)}, init ${(top.initiationProbability||0)*100|0}%)`,
                targetStartPos: top.start,
                regions: [
                    { name: '5UTR with uORF', type: '5UTR', start: fiveEnd, end: fiveStop },
                    { name: 'Renilla Luciferase', type: 'RLUC', start: fiveStop + 1, end: Math.min(fiveStop + 1000, sequence.length) },
                    { name: 'Linker', type: 'LINKER', start: Math.min(fiveStop + 1001, sequence.length), end: Math.min(fiveStop + 1030, sequence.length) },
                    { name: 'Firefly Luciferase', type: 'FLUC', start: Math.min(fiveStop + 1031, sequence.length), end: Math.min(fiveStop + 2950, sequence.length) }
                ]
            });
        }

        // Canonical CDS-based construct
        if (features.canonical && features.canonical.cds) {
            const cds = features.canonical.cds;
            suggestions.push({
                type: 'Main CDS Reporter',
                reason: `Targets canonical start at ${cds.start} (Kozak ${(features.canonical.start?.initiationProbability||0)*100|0}%)`,
                targetStartPos: cds.start,
                regions: [
                    { name: '5UTR', type: '5UTR', start: 1, end: Math.max(1, cds.start - 1) },
                    { name: 'Renilla Luciferase', type: 'RLUC', start: cds.start + 1, end: Math.min(cds.start + 1000, sequence.length) },
                    { name: 'Linker', type: 'LINKER', start: Math.min(cds.start + 1001, sequence.length), end: Math.min(cds.start + 1030, sequence.length) },
                    { name: 'Firefly Luciferase', type: 'FLUC', start: Math.min(cds.start + 1031, sequence.length), end: Math.min(cds.start + 2950, sequence.length) }
                ]
            });
        }
    } else {
        // Fallback heuristics if features missing
        const uorfs = startCodons.filter(sc => sc.pos < sequence.length / 3 && sc.orfLength < 300 && sc.orfLength > 9);
        if (uorfs.length > 0) {
            const lastUorf = Math.max(...uorfs.map(u => u.stopPos));
            suggestions.push({
                type: '5UTR',
                start: 1,
                end: Math.min(lastUorf + 50, Math.floor(sequence.length / 3)),
                reason: `Contains ${uorfs.length} uORF(s) - good for testing translational regulation`,
                regions: [
                    { name: '5UTR with uORFs', type: '5UTR', start: 1, end: lastUorf + 50 },
                    { name: 'Renilla Luciferase', type: 'RLUC', start: lastUorf + 51, end: lastUorf + 1050 },
                    { name: 'Linker', type: 'LINKER', start: lastUorf + 1051, end: lastUorf + 1080 },
                    { name: 'Firefly Luciferase', type: 'FLUC', start: lastUorf + 1081, end: Math.min(lastUorf + 2900, sequence.length) }
                ]
            });
        }
    }

    return suggestions;
}

// Export functions
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        calculateTreeLayout,
        drawTreeLayout,
        suggestConstructDesign
    };
} else {
    window.RDGViz = {
        calculateTreeLayout,
        drawTreeLayout,
        suggestConstructDesign
    };
}
