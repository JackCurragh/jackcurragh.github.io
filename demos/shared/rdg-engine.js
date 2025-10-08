/**
 * RDG Engine - Ribosome Decision Graph Calculation Module
 *
 * This module provides shared functionality for modeling translation initiation,
 * ribosome scanning, and alternative start site usage.
 *
 * Used by:
 * - RDG Demo (visualization and exploration)
 * - Western Blot Predictor (reporter construct analysis)
 */

// ============================================================================
// CONSTANTS
// ============================================================================

const RDG_CONSTANTS = {
    START_CODONS: {
        AUG: { isCanonical: true, name: 'AUG' },
        CUG: { isCanonical: false, name: 'CUG' },
        UUG: { isCanonical: false, name: 'UUG' },
        GUG: { isCanonical: false, name: 'GUG' },
        ACG: { isCanonical: false, name: 'ACG' },
        AUU: { isCanonical: false, name: 'AUU' },
        AUA: { isCanonical: false, name: 'AUA' },
        AUC: { isCanonical: false, name: 'AUC' }
    },
    STOP_CODONS: ['UAA', 'UAG', 'UGA'],

    // Kozak consensus positions relative to AUG at 0
    KOZAK_CONSENSUS: {
        '-3': 'G',  // Critical
        '+4': 'G'   // Critical
    }
};

// Default RDG model parameters
const DEFAULT_RDG_PARAMS = {
    baseP: 0.9,              // Base initiation probability for ideal AUG (90%)
    nearCognatePenalty: 0.5, // Multiplier for near-cognate starts
    distanceD0: 40,          // Distance decay constant (for ranking/auto-population only)
    gcBonus: 0.3,            // Bonus for high downstream GC (structure slows scanning)
    reinitiationBase: 0.3,   // Base reinitiation rate
    lengthL0: 100,           // uORF length decay constant
    spacingS0: 50            // Intercistronic spacing decay constant
};

// ============================================================================
// SEQUENCE ANALYSIS
// ============================================================================

/**
 * Find all start codons in a sequence
 * @param {string} sequence - RNA sequence (AUGC format)
 * @param {Array} readthroughStops - Array of stop positions that can read through (optional)
 * @returns {Array} Array of start codon objects with position, frame, type, etc.
 */
function findStartCodons(sequence, readthroughStops = []) {
    const startCodons = [];
    const startCodonList = Object.keys(RDG_CONSTANTS.START_CODONS);

    for (let i = 0; i < sequence.length - 2; i++) {
        const codon = sequence.substring(i, i + 3);

        if (startCodonList.includes(codon)) {
            const frame = i % 3;
            const isAUG = codon === 'AUG';

            // Find stop codon in this frame (respecting readthrough stops)
            let stopPos = sequence.length;
            for (let j = i + 3; j < sequence.length - 2; j += 3) {
                const stopCodon = sequence.substring(j, j + 3);
                if (RDG_CONSTANTS.STOP_CODONS.includes(stopCodon)) {
                    // Check if this stop can read through
                    const isReadthrough = readthroughStops.some(rt => rt.pos === j);
                    if (!isReadthrough) {
                        stopPos = j + 3; // Include the stop codon
                        break;
                    }
                }
            }

            // Calculate Kozak context score
            const kozakScore = calculateKozakScore(sequence, i);

            // Calculate downstream GC content
            const gcContent = calculateDownstreamGC(sequence, i);

            startCodons.push({
                pos: i,
                codon: codon,
                frame: frame,
                isAUG: isAUG,
                stopPos: stopPos,
                orfLength: stopPos - i,
                kozakScore: kozakScore,
                gcContent: gcContent
            });
        }
    }

    return startCodons;
}

/**
 * Calculate Kozak context score (0 to 1)
 * Based on positions -3 and +4 relative to start codon
 */
function calculateKozakScore(sequence, startPos) {
    let score = 0;

    // Check -3 position (most important)
    if (startPos >= 3) {
        const minus3 = sequence[startPos - 3];
        if (minus3 === 'G' || minus3 === 'A') score += 0.6; // G is optimal, A is good
    }

    // Check +4 position (also important)
    if (startPos + 3 < sequence.length) {
        const plus4 = sequence[startPos + 3];
        if (plus4 === 'G') score += 0.4; // G is optimal
    }

    return score;
}

/**
 * Calculate downstream GC content in a window after the start codon
 * High GC = stronger secondary structure = slower scanning = more initiation
 */
function calculateDownstreamGC(sequence, startPos, windowSize = 30) {
    // Look DOWNSTREAM from the start codon (after position + 3)
    const windowStart = startPos + 3;
    const windowEnd = Math.min(sequence.length, windowStart + windowSize);
    const window = sequence.substring(windowStart, windowEnd);

    let gcCount = 0;
    for (let i = 0; i < window.length; i++) {
        if (window[i] === 'G' || window[i] === 'C') gcCount++;
    }

    return window.length > 0 ? gcCount / window.length : 0.5;
}

/**
 * Find all stop codons in a sequence
 * @param {string} sequence - RNA sequence (AUGC format)
 * @returns {Array} Array of stop codon objects with position, frame, and codon type
 */
function findStopCodons(sequence) {
    const stopCodons = [];

    for (let i = 0; i < sequence.length - 2; i++) {
        const codon = sequence.substring(i, i + 3);
        if (RDG_CONSTANTS.STOP_CODONS.includes(codon)) {
            stopCodons.push({
                pos: i,
                codon: codon,
                frame: i % 3
            });
        }
    }

    return stopCodons;
}

/**
 * Identify canonical and predicted features within an RNA sequence
 * @param {string} sequence - RNA sequence (AUGC)
 * @param {Object} options - Additional options
 * @param {Array} options.readthroughStops - Optional readthrough stop annotations
 * @param {Object} options.params - Optional RDG parameter overrides
 * @returns {Object} Feature annotation object
 */
function identifyFeatures(sequence, options = {}) {
    const params = options.params || DEFAULT_RDG_PARAMS;
    const readthroughStops = options.readthroughStops || [];

    const startCodons = findStartCodons(sequence, readthroughStops).map((start, idx) => {
        const initiationProbability = calculateInitiationProbability(
            start.codon,
            start.kozakScore,
            start.gcContent,
            params
        );

        return {
            ...start,
            index: idx,
            initiationProbability
        };
    });

    const stopCodons = findStopCodons(sequence);
    const sortedStarts = [...startCodons].sort((a, b) => a.pos - b.pos);

    const canonicalStart = sortedStarts
        .filter(start => start.isAUG)
        .sort((a, b) => {
            if (b.initiationProbability !== a.initiationProbability) {
                return b.initiationProbability - a.initiationProbability;
            }
            return b.orfLength - a.orfLength;
        })[0] || sortedStarts[0] || null;

    const canonicalCDS = canonicalStart ? {
        start: canonicalStart.pos,
        end: canonicalStart.stopPos,
        frame: canonicalStart.frame,
        startCodon: canonicalStart.codon,
        initiationProbability: canonicalStart.initiationProbability
    } : null;

    const utr5 = canonicalStart ? {
        start: 0,
        end: canonicalStart.pos
    } : null;

    const utr3 = canonicalStart ? {
        start: canonicalStart.stopPos,
        end: sequence.length
    } : null;

    const uorfs = [];
    if (canonicalStart) {
        sortedStarts.forEach(start => {
            const isUpstream = start.pos < canonicalStart.pos;
            const terminatesBeforeCDS = start.stopPos <= canonicalStart.pos;

            if (isUpstream && terminatesBeforeCDS) {
                uorfs.push({
                    start: start.pos,
                    end: start.stopPos,
                    frame: start.frame,
                    startCodon: start.codon,
                    kozakScore: start.kozakScore,
                    initiationProbability: start.initiationProbability,
                    orfLength: start.orfLength
                });
            }
        });
    }

    const reinitiationSites = [];
    for (let i = 1; i < sortedStarts.length; i++) {
        const previous = sortedStarts[i - 1];
        const current = sortedStarts[i];
        const spacing = current.pos - previous.stopPos;

        if (spacing >= 0) {
            const reinitProb = calculateReinitiationProbability(previous.orfLength, spacing, params);
            if (reinitProb > 0.05) {
                reinitiationSites.push({
                    donorStart: previous.pos,
                    donorStop: previous.stopPos,
                    acceptorStart: current.pos,
                    spacing,
                    donorORFLength: previous.orfLength,
                    probability: reinitProb
                });
            }
        }
    }

    return {
        canonical: {
            start: canonicalStart || null,
            cds: canonicalCDS,
            utr5,
            utr3
        },
        predicted: {
            startCodons: sortedStarts,
            stopCodons,
            uorfs,
            reinitiationSites
        }
    };
}

/**
 * Find the next in-frame stop codon from a given position
 * @param {string} sequence - RNA sequence
 * @param {number} startPos - Starting position
 * @param {number} frame - Reading frame (0, 1, or 2)
 * @param {Array} readthroughStops - Stops that can read through
 * @returns {number} Position of next stop (or sequence length if none found)
 */
function findNextStopInFrame(sequence, startPos, frame, readthroughStops = []) {
    // Start searching from next codon in frame
    const searchStart = startPos + 3 - (startPos % 3) + frame;

    for (let j = searchStart; j < sequence.length - 2; j += 3) {
        if (j % 3 !== frame) continue; // Ensure we're in the right frame

        const stopCodon = sequence.substring(j, j + 3);
        if (RDG_CONSTANTS.STOP_CODONS.includes(stopCodon)) {
            // Check if this stop can read through
            const isReadthrough = readthroughStops.some(rt => rt.pos === j);
            if (!isReadthrough) {
                return j + 3; // Include the stop codon
            }
        }
    }

    return sequence.length;
}

/**
 * Apply a frameshift to a position and return new frame
 * @param {number} position - Current position
 * @param {number} shift - Frameshift amount (-1 or +1)
 * @returns {Object} New position and frame
 */
function applyFrameshift(position, shift) {
    const newPos = position + shift;
    const newFrame = newPos % 3;
    return { position: newPos, frame: newFrame };
}

// ============================================================================
// PROBABILITY CALCULATIONS
// ============================================================================

/**
 * Calculate initiation probability for a start codon
 * P_init = P₀ × f_kozak × f_codon × f_structure
 *
 * @param {string} startCodon - The start codon (AUG, CUG, etc.)
 * @param {number} kozakScore - Kozak context score (0-1)
 * @param {number} gcContent - Downstream GC content (0-1)
 * @param {Object} params - RDG model parameters (optional)
 * @returns {number} Initiation probability (0-1)
 */
function calculateInitiationProbability(startCodon, kozakScore, gcContent, params = DEFAULT_RDG_PARAMS) {
    // f_kozak: Kozak context factor (0.5 to 1.0)
    const fKozak = 0.5 + (kozakScore * 0.5);

    // f_codon: Start codon identity (1.0 for AUG, penalty for near-cognates)
    const isAUG = startCodon === 'AUG';
    const fCodon = isAUG ? 1.0 : params.nearCognatePenalty;

    // f_structure: Bonus for high downstream GC (structure slows scanning, increases initiation)
    // GC > 0.5 gives a bonus, GC < 0.5 is neutral
    const fStructure = 1.0 + (Math.max(0, gcContent - 0.5) * params.gcBonus);

    // Combined probability (NO DISTANCE FACTOR - distance only used for auto-population ranking)
    const pStart = params.baseP * fKozak * fCodon * fStructure;

    return Math.max(0.01, Math.min(0.99, pStart)); // Clamp between 1% and 99%
}

/**
 * Calculate reinitiation probability after completing a translon
 * Based on uORF length and spacing to next start codon
 *
 * From Gunisova et al. 2018, Hinnebusch 2014:
 * - Short uORFs (<35 codons) allow efficient reinitiation
 * - Long uORFs (>35 codons) suppress reinitiation
 * - Greater spacing increases reinitiation probability
 */
function calculateReinitiationProbability(uorfLength, spacing, params = DEFAULT_RDG_PARAMS) {
    // Length factor: shorter uORFs allow better reinitiation
    const fLength = Math.exp(-uorfLength / params.lengthL0);

    // Spacing factor: greater distance allows 60S rejoining
    const fSpacing = 1 - Math.exp(-spacing / params.spacingS0);

    // Combined reinitiation probability
    const pReinit = params.reinitiationBase * fLength * fSpacing;

    return Math.max(0.01, Math.min(0.95, pReinit));
}

/**
 * Calculate distance-based selection probability (for auto-population ranking only)
 * Favors start codons closer to 5' end
 */
function calculateDistanceWeight(position, params = DEFAULT_RDG_PARAMS) {
    return Math.exp(-position / params.distanceD0);
}

// ============================================================================
// PROTEIN CALCULATIONS
// ============================================================================

/**
 * Calculate molecular weight of a protein from ORF length
 * Average amino acid weight ≈ 110 Da
 *
 * @param {number} orfLengthNt - ORF length in nucleotides
 * @returns {number} Molecular weight in kDa
 */
function calculateProteinMW(orfLengthNt) {
    const aaCount = Math.floor(orfLengthNt / 3);
    const mwDa = aaCount * 110; // Average AA weight
    return mwDa / 1000; // Convert to kDa
}

// ============================================================================
// FLUX MODEL
// ============================================================================

/**
 * Compute scanning flux-based abundances for each start codon along the mRNA.
 * Simplified 1D model:
 *  - availableFlux starts at 1.0 at the 5' end
 *  - at start i with probability p_i, initiatedFlux = availableFlux * p_i
 *  - residualFlux to next start = availableFlux * (1 - p_i) + initiatedFlux * r_i
 *    where r_i is reinit probability based on ORF length and spacing to next start
 *  - abundance assigned to start i is initiatedFlux
 *
 * @param {Array} starts - start codon annotations (pos, stopPos, frame, etc.)
 * @param {Object} params - RDG params
 * @returns {Map<number, number>} Map from start.pos to abundance (0..1)
 */
function computeTranslonFlux(starts, params = DEFAULT_RDG_PARAMS) {
    const sorted = [...starts].sort((a, b) => a.pos - b.pos);
    const flux = new Map();

    let available = 1.0;
    for (let i = 0; i < sorted.length; i++) {
        const s = sorted[i];
        const p = s.initiationProbability !== undefined
            ? s.initiationProbability
            : calculateInitiationProbability(s.codon, s.kozakScore, s.gcContent, params);

        const initiated = available * Math.max(0, Math.min(1, p));

        flux.set(s.pos, initiated);

        const next = sorted[i + 1];
        const spacing = next ? Math.max(0, next.pos - s.stopPos) : 0;
        const r = calculateReinitiationProbability(s.orfLength, spacing, params);
        const reinitFlux = initiated * r;

        available = (available * (1 - p)) + reinitFlux;
        // Prevent drift outside [0,1]
        available = Math.max(0, Math.min(1, available));
    }

    // Normalize abundances to the maximum observed so visualizations are stable
    let max = 0;
    flux.forEach(v => { if (v > max) max = v; });
    if (max > 0) {
        const scaled = new Map();
        flux.forEach((v, k) => scaled.set(k, v / max));
        return scaled;
    }
    return flux;
}

/**
 * Translate nucleotide sequence to amino acid sequence
 * @param {string} ntSequence - Nucleotide sequence (must be multiple of 3)
 * @returns {string} Amino acid sequence (single letter codes)
 */
function translateSequence(ntSequence) {
    const codonTable = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    };

    let protein = '';
    for (let i = 0; i < ntSequence.length - 2; i += 3) {
        const codon = ntSequence.substring(i, i + 3);
        protein += codonTable[codon] || 'X'; // X for unknown
    }

    return protein;
}

// ============================================================================
// REPORTER-SPECIFIC CALCULATIONS
// ============================================================================

/**
 * Known reporter proteins with their molecular weights
 */
const REPORTER_PROTEINS = {
    RLUC: { name: 'Renilla Luciferase', mw: 36 }, // kDa
    FLUC: { name: 'Firefly Luciferase', mw: 61 },
    GFP: { name: 'Green Fluorescent Protein', mw: 27 },
    FLAG: { name: 'FLAG Tag', mw: 1 },
    HA: { name: 'HA Tag', mw: 1 }
};

/**
 * Reporter ORF library (synthetic coding sequences approximating target MWs)
 * Generates minimal AUG...(AA)...STOP sequences at desired length.
 */
const REPORTER_LIBRARY = (function() {
    function aaFromMw(mwKDa) {
        return Math.max(30, Math.round((mwKDa * 1000) / 110));
    }
    function synthesizeORF(aaLen) {
        // Start codon AUG, then varied safe codons (no internal stops, avoid near-cognate starts), then UAA stop
        // Avoid near-cognates: ACG, AUU, AUA, AUC, CUG, UUG, GUG
        const safeCodons = [
            'GCU','GCC','GCA','GCG',   // Ala
            'GAU','GAC',               // Asp
            'GAA','GAG',               // Glu
            'GGU','GGC','GGA','GGG',   // Gly
            'AAU','AAC',               // Asn (AAU/AAC safe; near-cognates are AUU/AUA/AUC)
            'AAA','AAG',               // Lys
            'ACU','ACC','ACA',         // Thr (exclude ACG near-cognate)
            'CAU','CAC',               // His
            'UAU','UAC',               // Tyr
            'UGU','UGC',               // Cys
            'CCU','CCC','CCA','CCG'    // Pro
        ];
        const bodyLen = Math.max(0, aaLen - 1);
        let body = '';
        for (let i = 0; i < bodyLen; i++) {
            body += safeCodons[i % safeCodons.length];
        }
        return 'AUG' + body + 'UAA';
    }
    function synthesizeORFNoStop(aaLen) {
        // Varied codons without terminal stop (avoid near-cognates)
        const safeCodons = [
            'GCU','GCC','GCA','GCG','GAU','GAC','GAA','GAG','GGU','GGC','GGA','GGG',
            'AAU','AAC','AAA','AAG','ACU','ACC','ACA','CAU','CAC','UAU','UAC','UGU','UGC','CCU','CCC','CCA','CCG'
        ];
        const bodyLen = Math.max(0, aaLen - 1);
        let body = '';
        for (let i = 0; i < bodyLen; i++) {
            body += safeCodons[i % safeCodons.length];
        }
        return 'AUG' + body; // no stop codon
    }
    function synthesizeORFWeakPlus4(aaLen) {
        // First codon after AUG starts with A (ACU) so +4 != G, rest alanines
        const after = 'ACU';
        const bodyLen = Math.max(0, aaLen - 2);
        let body = '';
        for (let i = 0; i < bodyLen; i++) body += 'GCU';
        return 'AUG' + after + body + 'UAA';
    }
    const rlucAA = aaFromMw(REPORTER_PROTEINS.RLUC.mw);
    const flucAA = aaFromMw(REPORTER_PROTEINS.FLUC.mw);
    const gfpAA = aaFromMw(REPORTER_PROTEINS.GFP.mw);
    // If the host page provides canonical reporter CDS in RNA (AUGC), prefer them.
    // window.REPORTER_CDS = { RLUC: 'AUG...', FLUC: 'AUG...' }
    try {
        if (typeof window !== 'undefined' && window.REPORTER_CDS) {
            const toRNA = (s) => (s||'').trim().toUpperCase().replace(/[^ACGTU]/g,'').replace(/T/g,'U');
            const cds = {
                RLUC: toRNA(window.REPORTER_CDS.RLUC),
                FLUC: toRNA(window.REPORTER_CDS.FLUC)
            };
            const hasRLUC = cds.RLUC && cds.RLUC.startsWith('AUG') && cds.RLUC.length >= 60;
            const hasFLUC = cds.FLUC && cds.FLUC.startsWith('AUG') && cds.FLUC.length >= 60;
            if (hasRLUC || hasFLUC) {
                const rluc = hasRLUC ? cds.RLUC : synthesizeORF(rlucAA);
                const fluc = hasFLUC ? cds.FLUC : synthesizeORF(flucAA);
                const rlucWeak = (function(){
                    // weaken +4 (first nt after AUG) if possible
                    if (rluc.length >= 4) {
                        const plus4 = rluc[3];
                        const repl = (plus4 === 'G') ? 'A' : plus4; // avoid strong G if present
                        return rluc.substring(0,3) + repl + rluc.substring(4);
                    }
                    return rluc;
                })();
                const rlucNoStop = (function(){
                    const tail = rluc.slice(-3);
                    if (['UAA','UAG','UGA'].includes(tail)) return rluc.slice(0,-3);
                    return rluc;
                })();
                return {
                    RLUC: rluc,
                    RLUC_NO_STOP: rlucNoStop,
                    RLUC_WEAK: rlucWeak,
                    FLUC: fluc,
                    GFP: synthesizeORF(gfpAA),
                    LINKER: (function(){
                        // 10x Gly/Ser flexible linker
                        let s='';
                        for(let i=0;i<10;i++) s += 'GGUGGUAGU'; // GGU=Gly, AGU=Ser (RNA)
                        return s;
                    })()
                };
            }
        }
    } catch {}

    return {
        RLUC: synthesizeORF(rlucAA),
        RLUC_NO_STOP: synthesizeORFNoStop(rlucAA),
        RLUC_WEAK: synthesizeORFWeakPlus4(rlucAA),
        FLUC: synthesizeORF(flucAA),
        GFP: synthesizeORF(gfpAA),
        LINKER: (function(){
            // 10x Gly/Ser flexible linker
            let s='';
            for(let i=0;i<10;i++) s += 'GGUGGUAGU'; // GGU=Gly, AGU=Ser (RNA)
            return s;
        })()
    };
})();

// Return reporter CDS dynamically, preferring canonical sequences if provided on window.REPORTER_CDS
function getReporterCDS(type) {
    try {
        const w = (typeof window !== 'undefined') ? window : null;
        const toRNA = (s) => (s||'').trim().toUpperCase().replace(/[^ACGTU]/g,'').replace(/T/g,'U');
        const stops = new Set(['UAA','UAG','UGA']);
        const deriveWeakPlus4 = (seq) => {
            if (!seq || seq.length < 4) return seq;
            const plus4 = seq[3];
            const repl = (plus4 === 'G') ? 'A' : plus4; // avoid strong G
            return seq.substring(0,3) + repl + seq.substring(4);
        };
        const stripTerminalStop = (seq) => {
            if (!seq || seq.length < 3) return seq;
            const tail = seq.slice(-3);
            if (stops.has(tail)) return seq.slice(0, -3);
            return seq;
        };
        if (w && w.REPORTER_CDS) {
            const rluc = toRNA(w.REPORTER_CDS.RLUC);
            const fluc = toRNA(w.REPORTER_CDS.FLUC);
            if ((type === 'RLUC' || type === 'RLUC_NO_STOP' || type === 'RLUC_WEAK') && rluc && rluc.startsWith('AUG')) {
                if (type === 'RLUC_WEAK') return deriveWeakPlus4(rluc);
                if (type === 'RLUC_NO_STOP') return stripTerminalStop(rluc);
                return rluc;
            }
            if (type === 'FLUC' && fluc && fluc.startsWith('AUG')) {
                return fluc;
            }
        }
    } catch (e) { /* ignore */ }

    // Fallback to static library
    return REPORTER_LIBRARY[type];
}

/**
 * Assemble a construct sequence from user regions.
 * - 5UTR/3UTR/CUSTOM slices are taken from baseSequence
 * - RLUC/FLUC/GFP use synthetic reporter ORFs from library
 * - LINKER uses a short flexible linker sequence
 * Returns assembled sequence and region coordinates in assembled space.
 */
function assembleConstruct(baseSequence, regions) {
    let assembled = '';
    const assembledRegions = [];

    function appendRegion(type, seq) {
        const start = assembled.length;
        assembled += seq;
        const end = assembled.length;
        assembledRegions.push({ type, start, end, assembled: true, color: (typeof window !== 'undefined' && window.REGION_COLORS && window.REGION_COLORS[type]) ? window.REGION_COLORS[type] : undefined });
    }

    for (const region of regions) {
        const t = region.type;
        if (t === '5UTR' || t === '3UTR' || t === 'CUSTOM') {
            const s = Math.max(0, (region.start || 1) - 1);
            const e = Math.min(baseSequence.length, region.end || baseSequence.length);
            if (e > s) appendRegion(t, baseSequence.substring(s, e));
        } else if (t === 'PAD_G3') {
            appendRegion(t, 'GGG');
        } else if (t === 'RLUC' || t === 'RLUC_NO_STOP' || t === 'RLUC_WEAK' || t === 'FLUC' || t === 'GFP') {
            if (t === 'GFP') {
                appendRegion(t, REPORTER_LIBRARY.GFP);
            } else {
                appendRegion(t, getReporterCDS(t));
            }
        } else if (t === 'LINKER') {
            appendRegion(t, REPORTER_LIBRARY.LINKER);
        } else {
            // Unknown: treat as slice
            const s = Math.max(0, (region.start || 1) - 1);
            const e = Math.min(baseSequence.length, region.end || baseSequence.length);
            if (e > s) appendRegion(t, baseSequence.substring(s, e));
        }
    }

    return { sequence: assembled, regions: assembledRegions };
}

function canonicalReporter(type) {
    if (type === 'RLUC' || type === 'RLUC_WEAK' || type === 'RLUC_NO_STOP') return 'RLUC';
    if (type === 'FLUC') return 'FLUC';
    if (type === 'GFP') return 'GFP';
    return null;
}

/**
 * Calculate expected protein products from a bicistronic reporter construct
 * @param {Array} translons - Array of translon objects (from RDG analysis)
 * @param {Object} constructMap - Map of ORF positions to reporter names
 * @returns {Array} Array of protein product objects with MW and abundance
 */
function predictProteinProducts(translons, constructMap) {
    const products = [];

    translons.forEach(translon => {
        // Determine which reporters this translon would produce
        const reporters = [];

        // Check overlap with reporter regions using frame-aware, lenient criteria
        for (const [orfType, region] of Object.entries(constructMap)) {
            const rep = canonicalReporter(orfType);
            if (!rep || !(rep in REPORTER_PROTEINS)) continue;

            const regionStart = region.start;
            const regionEnd = region.end;
            const regionLen = Math.max(0, regionEnd - regionStart);
            const regionFrame = regionStart % 3;

            // Overlap calculation
            const overlapStart = Math.max(translon.startNt, regionStart);
            const overlapEnd = Math.min(translon.endNt, regionEnd);
            const overlapLen = Math.max(0, overlapEnd - overlapStart);

            const startsInsideRegion = translon.startNt >= regionStart && translon.startNt < regionEnd;
            const spansMostOfRegion = overlapLen >= 0.5 * regionLen;
            const frameMatches = translon.frame === regionFrame;

            if (frameMatches && (startsInsideRegion || spansMostOfRegion)) {
                reporters.push(rep);
            }
        }

        // Calculate MW from ORF length for realism (handles fusions automatically)
        const totalMW = calculateProteinMW(translon.endNt - translon.startNt);

        const detectable = reporters.length > 0;
        products.push({
            name: translon.name,
            reporters: reporters,
            mw: totalMW,
            abundance: (translon.predictedAbundance ?? translon.predictedLUC ?? translon.probability) || 0,
            startPos: translon.startNt,
            endPos: translon.endNt,
            detectable
        });
    });

    return products;
}

/**
 * Build translon objects with RDG-derived abundance and protein size predictions
 * @param {string} sequence - RNA sequence (AUGC)
 * @param {Object} options - Additional options
 * @param {number} options.limit - Maximum number of translons to return
 * @param {Array} options.startCodons - Pre-computed start codon annotations
 * @param {Object} options.features - Pre-computed feature annotations
 * @param {Object} options.params - RDG parameter overrides
 * @returns {Array} Array of translon objects
 */
function buildTranslons(sequence, options = {}) {
    const limit = typeof options.limit === 'number' ? options.limit : 12;
    const params = options.params || DEFAULT_RDG_PARAMS;

    const featureData = options.features || identifyFeatures(sequence, { params });
    const canonicalStartPos = featureData.canonical.start ? featureData.canonical.start.pos : null;
    const uorfStartLookup = new Set(featureData.predicted.uorfs.map(u => u.start));

    const startAnnotations = options.startCodons && options.startCodons.length > 0
        ? options.startCodons.map((start, idx) => ({
            ...start,
            initiationProbability: start.initiationProbability || calculateInitiationProbability(
                start.codon,
                start.kozakScore,
                start.gcContent,
                params
            ),
            index: typeof start.index === 'number' ? start.index : idx
        }))
        : featureData.predicted.startCodons;

    // Compute scanning flux-based abundances along the mRNA
    const fluxByPos = computeTranslonFlux(startAnnotations, params);

    const makeTranslon = (start, idx) => {
        const initiationProbability = start.initiationProbability || calculateInitiationProbability(
            start.codon,
            start.kozakScore,
            start.gcContent,
            params
        );

        const proteinSize = parseFloat(calculateProteinMW(start.orfLength).toFixed(2));

        let classification = 'downstream';
        if (canonicalStartPos !== null && start.pos === canonicalStartPos) {
            classification = 'canonical';
        } else if (uorfStartLookup.has(start.pos)) {
            classification = 'uORF';
        } else if (canonicalStartPos !== null && start.pos < canonicalStartPos) {
            classification = 'upstream';
        }

        return {
            pathId: `T${idx + 1}`,
            name: `T${idx + 1}`,
            startNt: start.pos,
            endNt: start.stopPos,
            stopNt: start.stopPos,
            frame: start.frame,
            startCodon: start.codon,
            kozakScore: start.kozakScore,
            gcContent: start.gcContent,
            isCanonicalStart: classification === 'canonical',
            isAUG: start.isAUG,
            classification,
            predictedProteinSize: proteinSize,
            predictedProteinSizeKd: proteinSize,
            predictedAbundance: parseFloat((fluxByPos.get(start.pos) || initiationProbability).toFixed(4)),
            predictedLUC: parseFloat((fluxByPos.get(start.pos) || initiationProbability).toFixed(4)),
            probability: initiationProbability,
            orfLength: start.orfLength,
            initiationProbability,
            sourceIndex: start.index
        };
    };

    const translons = startAnnotations.map((start, idx) => makeTranslon(start, idx));

    translons.sort((a, b) => {
        if (b.predictedAbundance !== a.predictedAbundance) {
            return b.predictedAbundance - a.predictedAbundance;
        }
        return a.startNt - b.startNt;
    });

    let limited = translons.slice(0, Math.max(0, limit));

    // Ensure critical starts (e.g., reporter starts) are included if provided
    const ensureStarts = Array.isArray(options.ensureStarts) ? options.ensureStarts : [];
    ensureStarts.forEach((pos) => {
        if (!limited.some(t => t.startNt === pos)) {
            const start = startAnnotations.find(s => s.pos === pos);
            if (start) {
                const t = makeTranslon(start, translons.length + 1);
                limited.push(t);
            }
        }
    });

    // Re-label path IDs to reflect sorted order
    limited.forEach((translon, idx) => {
        const pathLabel = `T${idx + 1}`;
        translon.pathId = pathLabel;
        translon.name = pathLabel;
    });

    return limited;
}

// ============================================================================
// EXPORTS
// ============================================================================

// Export all functions and constants
if (typeof module !== 'undefined' && module.exports) {
    // Node.js/CommonJS
    module.exports = {
        RDG_CONSTANTS,
        DEFAULT_RDG_PARAMS,
        findStartCodons,
        findStopCodons,
        identifyFeatures,
        findNextStopInFrame,
        applyFrameshift,
        calculateKozakScore,
        calculateDownstreamGC,
        calculateInitiationProbability,
        calculateReinitiationProbability,
        calculateDistanceWeight,
        calculateProteinMW,
        translateSequence,
        computeTranslonFlux,
        assembleConstruct,
        REPORTER_PROTEINS,
        predictProteinProducts,
        buildTranslons
    };
} else {
    // Browser - attach to window
    window.RDGEngine = {
        RDG_CONSTANTS,
        DEFAULT_RDG_PARAMS,
        findStartCodons,
        findStopCodons,
        identifyFeatures,
        findNextStopInFrame,
        applyFrameshift,
        calculateKozakScore,
        calculateDownstreamGC,
        calculateInitiationProbability,
        calculateReinitiationProbability,
        calculateDistanceWeight,
        calculateProteinMW,
        translateSequence,
        computeTranslonFlux,
        assembleConstruct,
        REPORTER_PROTEINS,
        predictProteinProducts,
        buildTranslons
    };
}
