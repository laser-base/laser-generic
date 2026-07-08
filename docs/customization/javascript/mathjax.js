// MathJax configuration for pymdownx.arithmatex (generic mode) and Jupyter notebooks.
// Enables inline \(...\)/$...$ and display \[...\]/$$...$$ delimiters: \(...\) and \[...\]
// are what arithmatex emits for regular Markdown pages, while $...$ and $$...$$ are the raw
// delimiters nbconvert leaves in notebook cells (mkdocs-jupyter renders notebooks outside the
// Markdown/arithmatex pipeline, so its own math never gets converted to \(...\)/\[...\]).
// No ignoreHtmlClass/processHtmlClass restriction: MathJax's default skipHtmlTags already
// excludes script/style/pre/code, and a class-based allowlist can't reach notebook math anyway
// - it sits inside an unclassed <p> wrapper one level below the notebook's own output class, so
// the allowlist match doesn't survive the extra nesting.
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"], ["$", "$"]],
    displayMath: [["\\[", "\\]"], ["$$", "$$"]],
    processEscapes: true,
    processEnvironments: true
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
