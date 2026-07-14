#!/bin/bash
# run_all.sh — one command: compute results + build dashboard
# Usage:
#   ./run_all.sh                                    # Plan4RES toyDataset (default path)
#   ./run_all.sh /path/to/toyDataset RoGonDor       # explicit Plan4RES dataset
#   ./run_all.sh --omvs OMVS_inflows.csv            # real OMVS data
set -e
PY="${PYTHON:-/usr/bin/python3}"
cd "$(dirname "$0")"

if [ "$1" == "--omvs" ]; then
    echo ">>> Running MOSSHOOS on real OMVS data: $2"
    $PY compare_plan4res.py --inflow "$2" --n_iter 100 --n_mc 300
else
    DATASET="${1:-$HOME/plan4res_work/data/toyDataset}"
    RES="${2:-RoGonDor}"
    echo ">>> Running MOSSHOOS on Plan4RES inputs: $DATASET ($RES)"
    $PY compare_plan4res.py --dataset "$DATASET" --reservoir "$RES" --n_iter 100 --n_mc 300
fi

echo ">>> Building dashboard..."
$PY build_dashboard.py
echo ""
echo ">>> DONE. Open dashboard.html in your browser."
echo ">>> Results: results/comparison_results.json"
