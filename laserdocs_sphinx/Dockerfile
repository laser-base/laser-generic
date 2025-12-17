# Build docs for LASER-GENERIC
FROM sphinxdoc/sphinx-latexpdf:latest

# Install system pandoc and jupyter
RUN apt-get update && apt-get install -y pandoc && rm -rf /var/lib/apt/lists/*

# Install extra Sphinx extensions you use
# Pin myst-parser<4.0 because myst-parser 4.x introduced breaking changes in configuration and linkify behavior that this docs build is not yet compatible with.
RUN pip install --no-cache-dir \
      "myst-parser[linkify]<4.0" \
      sphinx-autodoc-typehints \
      sphinx-rtd-theme \
      sphinx-autoapi \
      jupyter \
      nbconvert \
      myst-nb

WORKDIR /docs/sphinx

ENTRYPOINT ["scripts/entrypoint.sh"]
CMD ["latexpdf"]
