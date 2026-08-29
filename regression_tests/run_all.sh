#!/usr/bin/env bash

set -euo pipefail

MADVORO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THUNDER_ROOT="${MADVORO_ROOT}/regression_tests/THUNDER"
exec "${THUNDER_ROOT}/run_all.sh" \
  --thunder-config "${MADVORO_ROOT}/regression_tests/config.json" "$@"
