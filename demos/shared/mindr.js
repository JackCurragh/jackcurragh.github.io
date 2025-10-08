/**
 * MINDR Control Generation Helpers
 *
 * Generates default control constructs for a selected test construct based on
 * RDG feature annotations. Uses simple codon substitutions compatible with the
 * demo's mutation model (in-place 3-nt codon replacements).
 */

(function(global){
  function mutateCodon(sequence, pos, newCodon) {
    return sequence.substring(0, pos) + newCodon + sequence.substring(pos + 3);
  }

  function generateControls(sequence, features, target) {
    // target: { name, startPos } — startPos is nt index of start codon
    const controls = [];

    if (!features || !features.predicted || !Number.isInteger(target.startPos)) {
      return controls;
    }

    const targetStart = target.startPos;
    const targetStartObj = (features.predicted.startCodons || []).find(s => s.pos === targetStart);
    const targetStop = targetStartObj ? targetStartObj.stopPos : (features.canonical && features.canonical.cds ? features.canonical.cds.end : targetStart + 90);

    // Test construct (baseline) - no mutations
    controls.push({
      label: `${target.name} – Test`,
      mutations: {}
    });

    // Initiation Negative Control: mutate target AUG to AAG (abolish initiation)
    controls.push({
      label: `${target.name} – Initiation(-) AAG`,
      mutations: { [targetStart]: 'AAG' }
    });

    // Termination Negative Control: introduce early stop codon ~10 codons into target ORF
    const earlyStopPos = Math.min(targetStop - 3, targetStart + 3 * 10);
    controls.push({
      label: `${target.name} – Early Stop`,
      mutations: { [earlyStopPos]: 'UAA' }
    });

    // Positive Control: remove upstream uORF starts (AUG->AAG) prior to target
    const upstreamStarts = (features.predicted.startCodons || []).filter(s => s.pos < targetStart && s.isAUG);
    const posMut = {};
    upstreamStarts.forEach(s => { posMut[s.pos] = 'AAG'; });
    if (Object.keys(posMut).length > 0) {
      controls.push({
        label: `${target.name} – Upstream(-)`,
        mutations: posMut
      });
    }

    return controls;
  }

  // Export
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = { generateControls };
  } else {
    global.MINDR = { generateControls };
  }
})(typeof window !== 'undefined' ? window : global);
