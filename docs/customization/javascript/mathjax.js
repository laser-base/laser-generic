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
    ignoreHtmlClass: ".*",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  // Guards against the MathJax CDN script failing to load (network block,
  // ad-blocker, offline): without this check, a failed load leaves
  // MathJax.startup undefined, throwing here and killing this subscription
  // for the rest of the session, including on later page navigations.
  if (!window.MathJax?.startup) return;
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
