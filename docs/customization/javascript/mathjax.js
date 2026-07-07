// MathJax configuration for pymdownx.arithmatex (generic mode).
// Enables inline \(...\) and display \[...\] delimiters that arithmatex emits,
// restricts typesetting to .arithmatex spans, and re-typesets after Material's
// instant navigation swaps page content (otherwise math fails to render on
// client-side page changes).
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
