window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\\[", "\\\\]"], ["$$", "$$"]],
    processEscapes: true,
    processEnvironments: true,
    tags: "ams"
  },
  options: {
    processHtmlClass: 'arithmatex'
    ignoreHtmlClass: "tex2jax_ignore",
  }
};

document$.subscribe(() => { 
  MathJax.typesetPromise()
})
