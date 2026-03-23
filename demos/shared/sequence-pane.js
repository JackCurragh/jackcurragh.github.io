/**
 * SequencePane - render RNA with simple annotations and actions
 */
(function(global){
  function wrapText(s, w){
    const lines=[]; for (let i=0;i<s.length;i+=w) lines.push(s.slice(i,i+w)); return lines.join('\n');
  }
  function render(container, rna, opts={}){
    const wrap = opts.wrap || 60;
    const starts = (opts.highlight && opts.highlight.starts) || [];
    const stops = (opts.highlight && opts.highlight.stops) || [];
    container.innerHTML = '';
    const bar = document.createElement('div');
    const copyBtn = document.createElement('button'); copyBtn.textContent = 'Copy'; copyBtn.className='secondary-btn';
    const dlBtn = document.createElement('button'); dlBtn.textContent='Download FASTA'; dlBtn.className='secondary-btn'; dlBtn.style.marginLeft='8px';
    const wrapSel = document.createElement('select'); wrapSel.innerHTML = '<option>60</option><option>100</option><option>120</option>'; wrapSel.value=String(wrap);
    bar.append('Wrap: ', wrapSel, ' ', copyBtn, ' ', dlBtn);
    bar.style.marginBottom = '6px';
    container.appendChild(bar);

    const pre = document.createElement('pre'); pre.style.whiteSpace='pre'; pre.style.userSelect='text'; pre.style.margin='0';
    function paint(){ pre.textContent = wrapText(rna||'', parseInt(wrapSel.value,10)); }
    paint();
    container.appendChild(pre);

    copyBtn.onclick = ()=> navigator.clipboard.writeText(rna||'');
    dlBtn.onclick = ()=>{
      const blob = new Blob([`>sequence\n${wrapText(rna||'', parseInt(wrapSel.value,10))}\n`], {type:'text/plain'});
      const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download='sequence.fasta'; a.click(); URL.revokeObjectURL(a.href);
    };
    wrapSel.onchange = paint;
  }

  global.SequencePane = { render };
})(typeof window !== 'undefined' ? window : globalThis);

