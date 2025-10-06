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

    for (let i = 0; i < sequence.length - 2; i += 3) {
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

        // Check if this translon overlaps with known reporters
        for (const [orfName, region] of Object.entries(constructMap)) {
            if (translon.startNt <= region.start && translon.endNt >= region.end) {
                reporters.push(orfName);
            }
        }

        // Calculate MW based on included reporters
        let totalMW = 0;
        reporters.forEach(rep => {
            if (REPORTER_PROTEINS[rep]) {
                totalMW += REPORTER_PROTEINS[rep].mw;
            }
        });

        // If no known reporters, calculate from ORF length
        if (totalMW === 0) {
            totalMW = calculateProteinMW(translon.endNt - translon.startNt);
        }

        products.push({
            name: translon.name,
            reporters: reporters,
            mw: totalMW,
            abundance: translon.probability, // Could be extended with flux data
            startPos: translon.startNt,
            endPos: translon.endNt
        });
    });

    return products;
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
        findNextStopInFrame,
        applyFrameshift,
        calculateKozakScore,
        calculateDownstreamGC,
        calculateInitiationProbability,
        calculateReinitiationProbability,
        calculateDistanceWeight,
        calculateProteinMW,
        translateSequence,
        REPORTER_PROTEINS,
        predictProteinProducts
    };
} else {
    // Browser - attach to window
    window.RDGEngine = {
        RDG_CONSTANTS,
        DEFAULT_RDG_PARAMS,
        findStartCodons,
        findStopCodons,
        findNextStopInFrame,
        applyFrameshift,
        calculateKozakScore,
        calculateDownstreamGC,
        calculateInitiationProbability,
        calculateReinitiationProbability,
        calculateDistanceWeight,
        calculateProteinMW,
        translateSequence,
        REPORTER_PROTEINS,
        predictProteinProducts
    };
}
