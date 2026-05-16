#! /bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"

BASE_DIR="${1:-$PROJECT_ROOT/RESULTS_CACTUS}"
MAX_JOBS="${2:-150}"
CONFIG_DIR="${3:-$PROJECT_ROOT/DRAMsim3/configs/cactus}"

trap "pkill -f sim_dramsim3" SIGINT

cd "$PROJECT_ROOT"
./scripts/run.sh -cfg "$CONFIG_DIR" -out "$BASE_DIR" -lim "$MAX_JOBS"
