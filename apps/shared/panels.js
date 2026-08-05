/**
 * Panels - three-state pane controller (expanded/dense/collapsed)
 */
(function(global){
  function wirePanel(root, opts={}){
    const def = opts.default || 'dense';
    let state = def;
    const btns = root.querySelectorAll('[data-pane-state]');
    function apply(){
      root.classList.remove('pane--expanded','pane--dense','pane--collapsed');
      root.classList.add(`pane--${state}`);
    }
    btns.forEach(b => b.addEventListener('click', ()=>{ state = b.getAttribute('data-pane-state'); apply(); if (api._cb) api._cb(state); }));
    apply();
    const api = { setState(s){ state=s; apply(); if (api._cb) api._cb(state); }, getState(){ return state; }, onChange(fn){ api._cb=fn; } };
    return api;
  }
  global.Panels = { wirePanel };
})(typeof window !== 'undefined' ? window : globalThis);

