# sphinx/scripts/convert_notebooks.py

import pathlib
import shutil

here = pathlib.Path(__file__).resolve().parent
root = here.parent.parent  # → /docs

nb_root = root / "docs" / "tutorials" / "notebooks"
target_dir = here.parent / "source" / "converted" / "tutorials"

target_dir.mkdir(parents=True, exist_ok=True)

for nb in nb_root.glob("*.ipynb"):
    target = target_dir / nb.name
    print(f"Copying {nb} → {target}")
    shutil.copy2(nb, target)
