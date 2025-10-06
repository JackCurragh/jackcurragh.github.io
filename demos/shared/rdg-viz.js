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

    // Draw construct regions if provided
    if (constructRegions && constructRegions.length > 0) {
        constructRegions.forEach(region => {
            const x1 = ntToPixel(region.start - 1);
            const x2 = ntToPixel(region.end);

            // Draw on top bar
            ctx.fillStyle = region.color;
            ctx.globalAlpha = 0.15;
            ctx.fillRect(x1, 20, x2 - x1, 40);
            ctx.globalAlpha = 1;

            ctx.strokeStyle = region.color;
            ctx.lineWidth = 2;
            ctx.strokeRect(x1, 20, x2 - x1, 40);

            ctx.fillStyle = region.color;
            ctx.font = 'bold 10px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText(region.type, (x1 + x2) / 2, 50);

            // If it's a reporter (RLUC or FLUC), draw the ORF on the appropriate frame track
            if (region.type === 'RLUC' || region.type === 'FLUC') {
                const frame = (region.start - 1) % 3; // Calculate reading frame
                const frameY = [80, 110, 140][frame];

                // Draw thicker bar on the frame track to show reporter ORF
                ctx.fillStyle = region.color;
                ctx.globalAlpha = 0.4;
                ctx.fillRect(x1, frameY - 10, x2 - x1, 20);
                ctx.globalAlpha = 1;

                ctx.strokeStyle = region.color;
                ctx.lineWidth = 3;
                ctx.strokeRect(x1, frameY - 10, x2 - x1, 20);

                // Add label with frame info
                ctx.fillStyle = region.color;
                ctx.font = 'bold 11px sans-serif';
                ctx.textAlign = 'center';
                const labelText = `${region.type} (Frame ${frame})`;
                ctx.fillText(labelText, (x1 + x2) / 2, frameY - 15);

                // Add start/stop markers for the reporter ORF
                ctx.strokeStyle = region.color;
                ctx.lineWidth = 2;
                ctx.setLineDash([]);

                // Start marker
                ctx.beginPath();
                ctx.moveTo(x1, frameY - 12);
                ctx.lineTo(x1, frameY + 12);
                ctx.stroke();

                // Stop marker
                ctx.beginPath();
                ctx.moveTo(x2, frameY - 12);
                ctx.lineTo(x2, frameY + 12);
                ctx.stroke();

                // Add arrow to show direction (5' to 3')
                ctx.fillStyle = region.color;
                ctx.beginPath();
                ctx.moveTo(x2 - 5, frameY);
                ctx.lineTo(x2, frameY - 4);
                ctx.lineTo(x2, frameY + 4);
                ctx.closePath();
                ctx.fill();
            }
        });
    }

    // Draw frame tracks
    const frameY = [80, 110, 140];
    const frameBarHeight = 20;
    const frameLabels = ['Frame 0', 'Frame 1', 'Frame 2'];

    ctx.font = 'bold 11px sans-serif';
    ctx.textAlign = 'right';

    frameY.forEach((y, frame) => {
        ctx.fillStyle = FRAME_COLORS[frame];
        ctx.fillText(frameLabels[frame], margin - 10, y + 6);

        ctx.fillStyle = FRAME_COLORS[frame];
        ctx.globalAlpha = 0.2;
        ctx.fillRect(margin, y - frameBarHeight/2, seqEndX - margin, frameBarHeight);
        ctx.globalAlpha = 1;

        ctx.strokeStyle = FRAME_COLORS[frame];
        ctx.lineWidth = 1;
        ctx.strokeRect(margin, y - frameBarHeight/2, seqEndX - margin, frameBarHeight);
    });

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
            ctx.strokeStyle = '#999';
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(x1, edge.y1);
            ctx.lineTo(x2, edge.y2);
            ctx.stroke();
        } else if (edge.type === 'vertical_branch') {
            ctx.strokeStyle = '#000';
            ctx.lineWidth = 3;
            ctx.beginPath();
            ctx.moveTo(x1, edge.y1);
            ctx.lineTo(x2, edge.y2);
            ctx.stroke();
        } else if (edge.type === 'translation') {
            ctx.strokeStyle = FRAME_COLORS[edge.translon.frame];
            ctx.lineWidth = 12;
            ctx.lineCap = 'butt';
            ctx.beginPath();
            ctx.moveTo(x1, edge.y1);
            ctx.lineTo(x2, edge.y2);
            ctx.stroke();

            // Determine which reporters this translon produces
            let reporters = [];
            if (constructRegions) {
                constructRegions.forEach(region => {
                    if ((region.type === 'RLUC' || region.type === 'FLUC') &&
                        edge.translon.startNt <= region.start - 1 &&
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
            const labelY = edge.y1 - 8;

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
                ctx.moveTo(stopX, edge.y1 - 6);
                ctx.lineTo(stopX, edge.y1 + 6);
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
                        const regY = [80, 110, 140][regFrame];

                        // Dashed connection line from translation bar to reporter
                        ctx.strokeStyle = region.color;
                        ctx.lineWidth = 1;
                        ctx.setLineDash([3, 3]);
                        ctx.globalAlpha = 0.4;

                        // Connect from middle of translation bar to reporter
                        const midX = (Math.max(x1, regX1) + Math.min(x2, regX2)) / 2;
                        ctx.beginPath();
                        ctx.moveTo(midX, edge.y1 + 6);
                        ctx.lineTo(midX, regY);
                        ctx.stroke();

                        ctx.setLineDash([]);
                        ctx.globalAlpha = 1;
                    }
                });
            }
        } else if (edge.type === 'readthrough') {
            // Readthrough segment - draw with dotted line to show it's conditional
            ctx.strokeStyle = FRAME_COLORS[edge.translon.frame];
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

            // Add readthrough probability label
            if (edge.probability !== undefined) {
                ctx.fillStyle = '#f59e0b';
                ctx.font = 'bold 10px sans-serif';
                ctx.textAlign = 'center';
                const midX = (x1 + x2) / 2;
                ctx.fillText(`RT ${(edge.probability * 100).toFixed(0)}%`, midX, edge.y1 - 5);
            }
        } else if (edge.type === 'reinitiation') {
            console.log(`Drawing reinitiation: x1=${x1}, y1=${edge.y1}, x2=${x2}, y2=${edge.y2}, label=${edge.label}`);
            ctx.strokeStyle = '#8b5cf6';
            ctx.lineWidth = 2;
            ctx.setLineDash([8, 4]);
            ctx.globalAlpha = 0.6;
            ctx.beginPath();
            ctx.moveTo(x1, edge.y1);
            ctx.lineTo(x2, edge.y2);
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
            ctx.arc(x, node.y, 5, 0, Math.PI * 2);
            ctx.fill();
        } else if (node.type === 'readthrough_decision') {
            // Special marker for readthrough decision points
            const x = ntToPixel(node.x);
            ctx.strokeStyle = '#f59e0b';
            ctx.fillStyle = '#fef3c7';
            ctx.lineWidth = 3;
            ctx.beginPath();
            ctx.arc(x, node.y, 7, 0, Math.PI * 2);
            ctx.fill();
            ctx.stroke();
        }
    });
}

// Suggest construct design based on sequence analysis
function suggestConstructDesign(sequence, startCodons) {
    const suggestions = [];

    // Find potential uORFs (short ORFs in 5' region)
    const uorfs = startCodons.filter(sc =>
        sc.pos < sequence.length / 3 && // In first third of sequence
        sc.orfLength < 300 && // Short ORF
        sc.orfLength > 9 // At least 3 amino acids
    );

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

    // Find main CDS (longest ORF)
    const mainCDS = startCodons.reduce((longest, current) =>
        current.orfLength > longest.orfLength ? current : longest
    , startCodons[0]);

    if (mainCDS && mainCDS.orfLength > 300) {
        suggestions.push({
            type: 'REPORTER_FUSION',
            start: mainCDS.pos,
            end: mainCDS.stopPos,
            reason: `Main CDS detected (${Math.floor(mainCDS.orfLength/3)} aa) - create N-terminal reporter fusion`,
            regions: [
                { name: '5UTR', type: '5UTR', start: 1, end: mainCDS.pos - 1 },
                { name: 'Renilla Luciferase', type: 'RLUC', start: mainCDS.pos, end: mainCDS.pos + 999 },
                { name: 'Native CDS', type: 'CUSTOM', start: mainCDS.pos + 1000, end: mainCDS.stopPos }
            ]
        });
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
