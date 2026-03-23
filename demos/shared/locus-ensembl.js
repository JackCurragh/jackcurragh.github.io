/**
 * Locus Ensembl Utilities
 *
 * Minimal helpers to fetch a gene with all transcripts/exons and the
 * genomic sequence covering the locus. Designed to be used from static
 * demos; no build step required.
 */

(function(global){
  const ENSEMBL_SERVER = global.ENSEMBL_SERVER || 'https://rest.ensembl.org';
  const SPECIES_NAME = global.SPECIES_NAME || 'homo_sapiens';

  function delay(ms){ return new Promise(r => setTimeout(r, ms)); }

  async function fetchWithRetry(url, opts={}, tries=3){
    let lastErr = null; let attempt = 0;
    while (attempt < tries){
      try {
        const res = await fetch(url, opts);
        if (!res.ok){
          const body = await res.text().catch(()=> '');
          // Retry on transient 5xx
          if (res.status >= 500 && res.status < 600 && attempt < tries-1){
            await delay(300 * Math.pow(2, attempt)); attempt++; continue;
          }
          throw new Error(`Ensembl error ${res.status}: ${url}\n${body}`);
        }
        return res;
      } catch (e){
        lastErr = e;
        if (attempt < tries-1){ await delay(300 * Math.pow(2, attempt)); attempt++; continue; }
        throw lastErr;
      }
    }
    throw lastErr || new Error('Unknown fetch error');
  }

  async function jsonFetch(url) {
    const res = await fetchWithRetry(url, { headers: { 'Accept': 'application/json' } });
    return res.json();
  }

  async function textFetch(url, acceptType = 'text/plain') {
    const res = await fetchWithRetry(url, { headers: { 'Accept': acceptType } });
    return res.text();
  }

  function isEnsemblId(s) {
    return /^(ENSG|ENST|ENSMUSG|ENSMUST)/i.test(s || '');
  }

  // Normalize a symbol or ID to a Gene object (with Transcript/Exon via expand=1)
  async function fetchGeneExpanded(identifier) {
    // Try symbol first, then ID fallback
    let gene;
    if (!isEnsemblId(identifier)) {
      try {
        gene = await jsonFetch(`${ENSEMBL_SERVER}/lookup/symbol/${SPECIES_NAME}/${encodeURIComponent(identifier)}?expand=1`);
      } catch (e) {
        // Fall through to ID lookup below
      }
    }
    if (!gene) {
      const obj = await jsonFetch(`${ENSEMBL_SERVER}/lookup/id/${encodeURIComponent(identifier)}?expand=1`);
      // If a transcript was provided, climb to its parent gene
      if (obj.object_type === 'Transcript' && obj.Parent) {
        gene = await jsonFetch(`${ENSEMBL_SERVER}/lookup/id/${encodeURIComponent(obj.Parent)}?expand=1`);
      } else if (obj.object_type === 'Gene') {
        gene = obj;
      } else {
        // Try Parent if present
        if (obj.Parent) {
          gene = await jsonFetch(`${ENSEMBL_SERVER}/lookup/id/${encodeURIComponent(obj.Parent)}?expand=1`);
        } else {
          throw new Error(`Identifier did not resolve to a Gene or Transcript: ${identifier}`);
        }
      }
    }

    if (!gene || !Array.isArray(gene.Transcript)) {
      throw new Error(`No transcripts found for ${identifier}`);
    }
    return gene;
  }

  // Fetch genomic sequence for the gene locus. Returns uppercase DNA (ACGT).
  async function fetchLocusSequence(gene) {
    const chr = gene.seq_region_name || gene.seq_region || gene.seq_region_id;
    const start = gene.start;
    const end = gene.end;
    if (!chr || !start || !end) throw new Error('Gene lacks region coordinates');
    const region = `${chr}:${start}..${end}`; // strand handled later via reverse complement if needed
    const fasta = await textFetch(`${ENSEMBL_SERVER}/sequence/region/${SPECIES_NAME}/${region}`, 'text/x-fasta');
    const dna = fasta.split('\n').filter(l => !l.startsWith('>')).join('').toUpperCase().replace(/[^ACGT]/g,'N');
    return { dna, chr, start, end, strand: gene.strand };
  }

  function revComp(dna) {
    const comp = { A:'T', T:'A', C:'G', G:'C', N:'N' };
    return dna.split('').reverse().map(b => comp[b] || 'N').join('');
  }

  // Extract exon DNA from the locus sequence, respecting gene strand
  function exonSequenceFromLocus(locus, exon) {
    const s = Math.max(1, exon.start) - locus.start; // 0-based offset in locus
    const e = Math.min(locus.end, exon.end) - locus.start; // inclusive end
    const slice = locus.dna.substring(s, e + 1);
    if (locus.strand === -1) {
      return revComp(slice);
    }
    return slice;
  }

  // Utility: convert DNA to RNA (U) upper-case
  function dnaToRNA(s) { return (s||'').toUpperCase().replace(/T/g,'U'); }

  global.LocusEnsembl = {
    fetchGeneExpanded,
    fetchLocusSequence,
    exonSequenceFromLocus,
    dnaToRNA
  };
})(typeof window !== 'undefined' ? window : globalThis);
