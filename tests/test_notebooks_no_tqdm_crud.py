"""Gate: committed tutorial notebooks must not carry tqdm progress-bar noise.

`docs/tqdm_strip.py` strips these lines automatically inside
`docs/execute_notebooks.py` (so the CI executed_nbs artifact ships clean), but
contributors who execute notebooks locally and then commit the result can
still slip crud into `docs/tutorials/notebooks/`. This test fails the CI
matrix (`check` job in `github-actions.yml`) if any committed notebook has
tqdm noise in its outputs, prompting the author to run:

    python docs/tqdm_strip.py docs/tutorials/notebooks/*.ipynb

before re-committing.
"""

import sys
from pathlib import Path

import pytest

# Import from docs/ — same shape as other laser-generic test files that reach
# into the repo root.
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "docs"))
from tqdm_strip import count_tqdm_lines  # noqa: E402


NOTEBOOK_DIR = REPO_ROOT / "docs" / "tutorials" / "notebooks"


@pytest.mark.parametrize("notebook", sorted(NOTEBOOK_DIR.glob("*.ipynb")), ids=lambda p: p.name)
def test_notebook_has_no_tqdm_progress_lines(notebook):
    n = count_tqdm_lines(notebook)
    assert n == 0, (
        f"{notebook.name} has {n} tqdm-progress line(s) in its output cells. "
        "Run `python docs/tqdm_strip.py docs/tutorials/notebooks/*.ipynb` to strip, "
        "then re-commit."
    )
