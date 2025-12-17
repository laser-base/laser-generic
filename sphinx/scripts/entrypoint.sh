#!/bin/sh
set -e
python3 scripts/convert_docs.py
python3 scripts/convert_notebooks.py
python3 scripts/update_index_rst.py
exec make "$@"
