"""
harvest_test_docstrings.py
--------------------------
Collects docstrings from all test_* functions and unittest.TestCase methods
inside the tests/ package and generates a feature_test_report.md file at
the project root.

Usage:
    python harvest_test_docstrings.py
"""

import importlib
import inspect
import pkgutil
import sys
import textwrap
import unittest
from pathlib import Path


def harvest_docstrings(package="tests"):
    docs = []

    # Ensure parent directory (project root) is on sys.path so imports like "tests.foo" work
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    package_path = project_root / package
    for _, modname, _ in pkgutil.walk_packages([str(package_path)]):
        try:
            module = importlib.import_module(f"{package}.{modname}")
        except Exception as e:
            print(f"⚠️ Skipping {modname}: {e}")
            continue

        # --- collect module-level test functions ---
        for name, obj in inspect.getmembers(module, inspect.isfunction):
            if name.startswith("test_"):
                doc = inspect.getdoc(obj)
                if doc:
                    docs.append({"test": f"{modname}.{name}", "doc": textwrap.dedent(doc)})

        # --- collect test methods inside classes ---
        for cname, cls in inspect.getmembers(module, inspect.isclass):
            if issubclass(cls, unittest.TestCase):
                for name, method in inspect.getmembers(cls, inspect.isfunction):
                    if name.startswith("test_"):
                        doc = inspect.getdoc(method)
                        if doc:
                            docs.append(
                                {
                                    "test": f"{modname}.{cname}.{name}",
                                    "doc": textwrap.dedent(doc),
                                }
                            )
    return docs


if __name__ == "__main__":
    docs = harvest_docstrings("tests")

    out = "\n\n".join(f"### {d['test']}\n{d['doc']}" for d in docs)

    # Write to project root
    output_path = Path(__file__).resolve().parent.parent / "feature_test_report.md"
    output_path.write_text(out, encoding="utf-8")

    print(f"Wrote {output_path.relative_to(Path.cwd())}")
