# generate_docs_md.py
from pathlib import Path

source_dir = Path("src/laser/generic")
output_dir = Path("docs/reference")
exclude = {"SI.py", "SIR.py", "SIS.py", "SEIR.py", "SEIRS.py", "SIRS.py"}

output_dir.mkdir(parents=True, exist_ok=True)

for py in source_dir.glob("*.py"):
    if py.name in exclude or py.name == "__init__.py":
        continue
    modname = py.stem
    with open(output_dir / f"{modname}.md", "w") as f:
        f.write(f"# `{modname}.py`\n\n::: laser.generic.{modname}\n")
