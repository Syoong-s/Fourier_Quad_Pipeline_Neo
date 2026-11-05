#!/usr/bin/env bash
set -euo pipefail
FQ_RUNNER_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$FQ_RUNNER_DIR/common.sh"
fq_load "${1:?Usage: bash build.sh /path/cpppipeline.env}"
exec bash "$FQ_RUNNER_DIR/runtime.sh" "$FQ_ENV_FILE" build
