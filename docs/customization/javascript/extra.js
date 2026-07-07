// Activate the "Was this page helpful?" feedback widget.
//
// Material renders the feedback form with a `hidden` attribute and only ships
// the JS that reveals it as part of its analytics-provider integration
// (partials/integrations/analytics/google.html). We don't set
// extra.analytics.provider, so that JS never loads and the form stays hidden.
// This block reproduces Material's native behavior (un-hide, disable on submit,
// reveal the matching thank-you note) without pulling in Google Analytics.
window.document$.subscribe(function() {
    const feedback = document.forms.feedback;
    if (feedback === undefined) return;

    // Reveal the form (Material leaves it hidden by default)
    feedback.hidden = false;

    for (const button of feedback.querySelectorAll("[type=submit]")) {
        button.addEventListener("click", function(event) {
            event.preventDefault();

            const page = document.location.pathname;
            const value = this.getAttribute("data-md-value");

            // Record the rating. The Axiom/do11y wiring is added with the rest
            // of the analytics work; until then, send a CustomEvent so the
            // handler can be attached without touching this file.
            document.dispatchEvent(new CustomEvent("page-feedback", {
                detail: { page: page, value: value }
            }));

            // Disable further input and show the matching thank-you note
            feedback.firstElementChild.disabled = true;
            const note = feedback.querySelector(
                ".md-feedback__note [data-md-value='" + value + "']"
            );
            if (note) note.hidden = false;
        });
    }
})
