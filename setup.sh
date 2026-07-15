#!/usr/bin/env bash
# =============================================================================
# setup.sh - MOSSHOOS x Plan4RES
#
# Installs and verifies everything needed to run MOSSHOOS, with or without
# Plan4RES, with or without a commercial solver licence.
#
#   ./setup.sh              # detect, install deps, run all unit tests
#   ./setup.sh --check      # detect only, change nothing
#   ./setup.sh --no-tests   # install, skip the test suites
#
# MOSSHOOS itself needs only numpy/pandas/scipy/scikit-learn. Plan4RES is
# OPTIONAL: without it you still get the full MOSSHOOS pipeline and the
# reimplemented Plan4RES deterministic method (M-D) for comparison.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

RUN_TESTS=1
CHECK_ONLY=0
for a in "$@"; do
  case "$a" in
    --check)    CHECK_ONLY=1 ;;
    --no-tests) RUN_TESTS=0 ;;
    -h|--help)  sed -n '2,16p' "$0"; exit 0 ;;
  esac
done

ok(){   printf '  [ ok ] %s\n' "$1"; }
warn(){ printf '  [warn] %s\n' "$1"; }
bad(){  printf '  [FAIL] %s\n' "$1"; }
hdr(){  printf '\n== %s ==\n' "$1"; }

FAILED=0

# -----------------------------------------------------------------------------
hdr "1. Python"
# -----------------------------------------------------------------------------
# Pick an interpreter that actually has _ctypes (some source builds lack it,
# which breaks pandas with a confusing ModuleNotFoundError).
PY=""
for cand in "${PYTHON:-}" /usr/bin/python3 python3 python3.12 python3.11 python3.10; do
  [ -z "$cand" ] && continue
  command -v "$cand" >/dev/null 2>&1 || continue
  if "$cand" -c "import ctypes" >/dev/null 2>&1; then PY="$(command -v "$cand")"; break; fi
done
if [ -z "$PY" ]; then
  bad "no python3 with a working _ctypes module found"
  echo "       (a source-built python without libffi will fail on 'import pandas')"
  exit 1
fi
ok "$PY  ($("$PY" --version 2>&1))"
[ "$PY" != "$(command -v python3)" ] && \
  warn "default python3 is $(command -v python3) - use PYTHON=$PY or call $PY directly"

# -----------------------------------------------------------------------------
hdr "2. Python packages"
# -----------------------------------------------------------------------------
NEED=(numpy pandas scipy sklearn)
MISSING=()
for m in "${NEED[@]}"; do
  if "$PY" -c "import $m" >/dev/null 2>&1; then
    v=$("$PY" -c "import $m;print(getattr($m,'__version__','?'))" 2>/dev/null)
    ok "$m $v"
  else
    MISSING+=("$m"); bad "$m missing"
  fi
done
if [ ${#MISSING[@]} -gt 0 ] && [ "$CHECK_ONLY" -eq 0 ]; then
  echo "  installing: ${MISSING[*]}"
  PKGS=${MISSING[*]/sklearn/scikit-learn}
  # PEP 668 environments need the override; try the polite way first
  "$PY" -m pip install --user $PKGS >/dev/null 2>&1 || \
  "$PY" -m pip install --user --break-system-packages $PKGS >/dev/null 2>&1 || \
  { bad "pip install failed - install manually: $PY -m pip install --break-system-packages $PKGS"; FAILED=1; }
  for m in "${MISSING[@]}"; do
    "$PY" -c "import $m" >/dev/null 2>&1 && ok "$m now available" || { bad "$m still missing"; FAILED=1; }
  done
fi

# -----------------------------------------------------------------------------
hdr "3. Data"
# -----------------------------------------------------------------------------
if [ -f data/OMVS_inflows.csv ]; then
  n=$(($(wc -l < data/OMVS_inflows.csv) - 1))
  ok "data/OMVS_inflows.csv ($n months)"
else
  warn "data/OMVS_inflows.csv not found - only Plan4RES-input runs will work"
fi

# -----------------------------------------------------------------------------
hdr "4. Plan4RES (optional)"
# -----------------------------------------------------------------------------
P4R_OK=0
P4R_DIR="${PLAN4RES_DIR:-$HOME/plan4res/p4r-env}"
DS="${PLAN4RES_DATASET:-$HOME/plan4res_work/data/toyDataset}"
if [ -d "$P4R_DIR" ]; then
  ok "Plan4RES environment: $P4R_DIR"
  if [ -f "$DS/csv_simul/SS_SeasonalStorage.csv" ]; then
    ok "dataset: $DS"; P4R_OK=1
  else
    warn "no dataset at $DS (set PLAN4RES_DATASET=...)"
  fi
else
  warn "Plan4RES not installed - MOSSHOOS runs fine without it"
  echo "         (the deterministic Plan4RES method M-D is reimplemented in"
  echo "          src/sddp_clean.py, so the comparison still works)"
fi

# -----------------------------------------------------------------------------
hdr "5. Solvers"
# -----------------------------------------------------------------------------
# MOSSHOOS needs NO external solver - the SDDP and HJB are self-contained.
ok "MOSSHOOS: no external solver required (self-contained SDDP + HJB)"
if [ "$P4R_OK" -eq 1 ]; then
  U="$P4R_DIR/scripts/add-ons/.build/sms++/Umbrella"
  FOUND=""
  for s in HiGHS CPX SCIP GRB; do
    if find "$U" -name "${s}MILPSolver.cpp.o" -print -quit 2>/dev/null | grep -q . ; then
      FOUND="$FOUND $s"
    fi
  done
  [ -n "$FOUND" ] && ok "Plan4RES MILP solvers compiled:$FOUND" || warn "could not detect compiled MILP solvers"
  if [ -n "${GUROBI_HOME:-}" ] || command -v gurobi_cl >/dev/null 2>&1; then
    ok "Gurobi detected - Plan4RES SDDP should run"
  else
    warn "no Gurobi licence detected"
    echo "         Plan4RES's own SDDP currently FAILS on HiGHS-only builds:"
    echo "         'oneStepBackward: no linearization is available'"
    echo "         See plan4res_fixes/PLAN4RES_SDDP_BUG.md for the full diagnosis."
    echo "         This does NOT affect MOSSHOOS: we reimplement their"
    echo "         deterministic method (M-D) so the comparison runs regardless."
  fi
fi

# -----------------------------------------------------------------------------
hdr "6. Unit tests"
# -----------------------------------------------------------------------------
if [ "$RUN_TESTS" -eq 1 ] && [ "$CHECK_ONLY" -eq 0 ]; then
  for t in sddp_clean sde_convergence hjb_clean wavelet_galerkin; do
    printf '  %-20s ' "$t"
    if out=$("$PY" "src/$t.py" 2>&1); then
      if grep -q "ALL TESTS PASSED" <<< "$out"; then echo "PASS"
      else echo "FAIL"; grep -E "TEST .*(PASS|FAIL)" <<< "$out" | sed 's/^/        /'; FAILED=1; fi
    else
      echo "ERROR"; tail -3 <<< "$out" | sed 's/^/        /'; FAILED=1
    fi
  done
else
  warn "tests skipped"
fi

# -----------------------------------------------------------------------------
hdr "Result"
# -----------------------------------------------------------------------------
if [ "$FAILED" -eq 0 ]; then
  echo "  Ready."
  echo ""
  echo "  Run the pipeline:"
  echo "    PYTHON=$PY ./run_all.sh --omvs data/OMVS_inflows.csv"
  [ "$P4R_OK" -eq 1 ] && echo "    PYTHON=$PY ./run_all.sh $DS RoGonDor"
  echo "    $PY storage_sweep.py --inflow data/OMVS_inflows.csv"
  echo ""
  echo "  See the results:"
  echo "    $PY dashboard_server.py --port 8050     # live, re-run from the browser"
  echo "    $PY make_static.py                      # standalone dashboard.html"
  exit 0
else
  echo "  Some checks failed - see above."
  exit 1
fi
