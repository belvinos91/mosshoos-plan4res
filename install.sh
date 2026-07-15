#!/usr/bin/env bash
# =============================================================================
# install.sh - MOSSHOOS x Plan4RES : one-paste install and run (macOS / Linux)
#
#   bash install.sh
#
# Works from inside a clone, or standalone (it will clone for you).
# Creates an isolated venv so nothing touches your system Python.
# Ends by opening the dashboard in your browser.
# =============================================================================
set -uo pipefail

REPO="${MOSSHOOS_REPO:-https://github.com/belvinos91/mosshoos-plan4res.git}"
DIR="${MOSSHOOS_DIR:-mosshoos-plan4res}"
SKIP_RUN=0
for a in "$@"; do [ "$a" = "--no-run" ] && SKIP_RUN=1; done

say(){ printf '\n\033[1;34m==>\033[0m %s\n' "$1"; }
ok(){  printf '    \033[32mok\033[0m   %s\n' "$1"; }
bad(){ printf '    \033[31mFAIL\033[0m %s\n' "$1"; }

# --- 1. get the code ---------------------------------------------------------
if [ -f "compare_plan4res.py" ] && [ -d "src" ]; then
  say "Already inside the repository"
else
  say "Cloning $REPO"
  if [ -d "$DIR" ]; then
    ok "$DIR exists - updating"; cd "$DIR" && git pull --ff-only 2>/dev/null || true
  else
    command -v git >/dev/null 2>&1 || { bad "git not found. Install: https://git-scm.com/downloads"; exit 1; }
    if ! git clone "$REPO" "$DIR"; then
      bad "clone failed"
      echo "    This repository is private. You need either:"
      echo "      - to be added as a collaborator, then sign in when git asks, or"
      echo "      - a zip of the repository from the author."
      exit 1
    fi
    cd "$DIR"
  fi
fi
ok "$(pwd)"

# --- 2. find a usable python -------------------------------------------------
say "Python"
PY=""
for c in python3.12 python3.11 python3.10 python3 python; do
  command -v "$c" >/dev/null 2>&1 || continue
  # must be >=3.9 AND have a working _ctypes (some source builds lack libffi,
  # which fails later with a confusing 'No module named _ctypes' from pandas)
  "$c" -c "import sys,ctypes; sys.exit(0 if sys.version_info>=(3,9) else 1)" >/dev/null 2>&1 \
    && { PY="$(command -v "$c")"; break; }
done
if [ -z "$PY" ]; then
  bad "no usable Python 3.9+ found"
  echo "    macOS : brew install python3      (or https://www.python.org/downloads/)"
  echo "    Linux : sudo apt install python3 python3-venv python3-pip"
  exit 1
fi
ok "$PY  ($("$PY" --version 2>&1))"

# --- 3. isolated environment -------------------------------------------------
say "Creating an isolated environment (.venv)"
if [ ! -d .venv ]; then
  "$PY" -m venv .venv 2>/dev/null || {
    bad "venv creation failed - on Debian/Ubuntu: sudo apt install python3-venv"; exit 1; }
fi
# shellcheck disable=SC1091
source .venv/bin/activate
ok "activated (nothing outside .venv is modified)"

say "Installing dependencies"
python -m pip install --upgrade pip >/dev/null 2>&1
if python -m pip install -r requirements.txt >/tmp/mosshoos_pip.log 2>&1; then
  ok "numpy pandas scipy scikit-learn matplotlib"
else
  bad "pip install failed - see /tmp/mosshoos_pip.log"; tail -5 /tmp/mosshoos_pip.log; exit 1
fi

# --- 4. verify ---------------------------------------------------------------
say "Verifying the method (unit tests)"
FAIL=0
for t in sddp_clean sde_convergence hjb_clean wavelet_galerkin; do
  printf '    %-18s ' "$t"
  if out=$(python "src/$t.py" 2>&1) && grep -q "ALL TESTS PASSED" <<< "$out"; then
    printf '\033[32mPASS\033[0m\n'
  else
    printf '\033[31mFAIL\033[0m\n'; FAIL=1
    grep -E "TEST .*(PASS|FAIL)" <<< "$out" | sed 's/^/           /'
  fi
done
[ "$FAIL" -eq 0 ] && ok "all components validated" || bad "some tests failed (see above)"

# --- 5. run ------------------------------------------------------------------
if [ "$SKIP_RUN" -eq 0 ]; then
  say "Running the pipeline on real OMVS data (about a minute)"
  python compare_plan4res.py --inflow data/OMVS_inflows.csv \
        --n_iter 150 --n_mc 300 > /tmp/mosshoos_run.log 2>&1 \
    && ok "pipeline complete" || { bad "pipeline failed"; tail -8 /tmp/mosshoos_run.log; }
  grep -E "^  (MS|MS0|MD|MD0) " /tmp/mosshoos_run.log | sed 's/^/    /'

  if [ ! -f results/storage_sweep.json ]; then
    say "Storage-sensitivity sweep (a few minutes)"
    python storage_sweep.py --inflow data/OMVS_inflows.csv \
          --out results/storage_sweep.json > /tmp/mosshoos_sweep.log 2>&1 \
      && ok "sweep complete" || bad "sweep failed - see /tmp/mosshoos_sweep.log"
  fi

  say "Building the dashboard"
  python make_static.py 2>&1 | sed 's/^/    /'
  D="$(pwd)/dashboard.html"
  if [ -f "$D" ]; then
    command -v open    >/dev/null 2>&1 && open "$D"
    command -v xdg-open >/dev/null 2>&1 && xdg-open "$D" >/dev/null 2>&1
    ok "opened $D"
  fi
fi

# --- done --------------------------------------------------------------------
cat <<EOF

  ------------------------------------------------------------------
  Done.  Everything below assumes you are in $(pwd)
  and have run:   source .venv/bin/activate

    open dashboard.html                       standalone, offline, shareable
    python dashboard_server.py --port 8050    live: change params, re-run
                                              then open http://localhost:8050

    python compare_plan4res.py --inflow data/OMVS_inflows.csv --n_iter 300
    python storage_sweep.py    --inflow data/OMVS_inflows.csv
    python src/sddp_clean.py                  re-verify any component

  Results: results/comparison_results.json  -  RESULTS.md
  ------------------------------------------------------------------

EOF
