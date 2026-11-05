#!/usr/bin/env bash
# Shared configuration loader. Only source env files you trust.
set -euo pipefail

fq_load() {
    FQ_ENV_FILE="$(cd -- "$(dirname -- "$1")" && pwd)/$(basename -- "$1")"
    # shellcheck disable=SC1090
    source "$FQ_ENV_FILE"
}

fq_check() {
    python3 "$FQ_RUNNER_DIR/../fqctl.py" check --json "$FQ_JSON"
}
