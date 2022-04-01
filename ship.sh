#!/bin/bash

set -euo pipefail
tox
python3 -m build
python3 -m twine upload dist/*
