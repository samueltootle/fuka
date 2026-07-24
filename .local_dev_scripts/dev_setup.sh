#!/bin/bash

python3 -m venv .venv
source .venv/bin/activate

pip install --upgrade pip
pip install -r requirements.txt

git config core.hooksPath $HOME_KADATH/.githooks
chmod +x .githooks/pre-commit