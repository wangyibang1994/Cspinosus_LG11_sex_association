#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
python "$ROOT/workflow/reproduce_key_numbers.py" "$ROOT"
python "$ROOT/workflow/validate_release.py" "$ROOT"
