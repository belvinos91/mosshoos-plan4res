# =============================================================================
# install.ps1 - MOSSHOOS x Plan4RES : one-paste install and run (Windows)
#
#   powershell -ExecutionPolicy Bypass -File install.ps1
#   (or in the VS Code PowerShell terminal:  .\install.ps1 )
#
# Works from inside a clone, or standalone (it will clone for you).
# Creates an isolated .venv so nothing touches your system Python.
# Ends by opening the dashboard in your browser.
# =============================================================================
param([switch]$NoRun)

$ErrorActionPreference = "Continue"
$Repo = if ($env:MOSSHOOS_REPO) { $env:MOSSHOOS_REPO } else { "https://github.com/belvinos91/mosshoos-plan4res.git" }
$Dir  = if ($env:MOSSHOOS_DIR)  { $env:MOSSHOOS_DIR }  else { "mosshoos-plan4res" }

function Say($m) { Write-Host "`n==> $m" -ForegroundColor Cyan }
function Ok($m)  { Write-Host "    ok   $m" -ForegroundColor Green }
function Bad($m) { Write-Host "    FAIL $m" -ForegroundColor Red }

# --- 1. get the code ---------------------------------------------------------
if ((Test-Path "compare_plan4res.py") -and (Test-Path "src")) {
    Say "Already inside the repository"
} else {
    Say "Cloning $Repo"
    if (Test-Path $Dir) {
        Ok "$Dir exists - updating"
        Set-Location $Dir
        git pull --ff-only 2>$null
    } else {
        if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
            Bad "git not found. Install: https://git-scm.com/download/win"
            exit 1
        }
        git clone $Repo $Dir
        if ($LASTEXITCODE -ne 0) {
            Bad "clone failed"
            Write-Host "    This repository is private. You need either:"
            Write-Host "      - to be added as a collaborator, then sign in when git asks, or"
            Write-Host "      - a zip of the repository from the author."
            exit 1
        }
        Set-Location $Dir
    }
}
Ok (Get-Location).Path

# --- 2. find a usable python -------------------------------------------------
Say "Python"
$Py = $null
foreach ($c in @("python", "python3", "py")) {
    $cmd = Get-Command $c -ErrorAction SilentlyContinue
    if (-not $cmd) { continue }
    & $c -c "import sys,ctypes; sys.exit(0 if sys.version_info>=(3,9) else 1)" 2>$null
    if ($LASTEXITCODE -eq 0) { $Py = $c; break }
}
if (-not $Py) {
    Bad "no usable Python 3.9+ found"
    Write-Host "    Install from https://www.python.org/downloads/windows/"
    Write-Host "    IMPORTANT: tick 'Add Python to PATH' during installation."
    exit 1
}
Ok "$Py  ($(& $Py --version 2>&1))"

# --- 3. isolated environment -------------------------------------------------
Say "Creating an isolated environment (.venv)"
if (-not (Test-Path ".venv")) {
    & $Py -m venv .venv
    if ($LASTEXITCODE -ne 0) { Bad "venv creation failed"; exit 1 }
}
$VPy = Join-Path (Get-Location) ".venv\Scripts\python.exe"
if (-not (Test-Path $VPy)) { Bad "venv python not found at $VPy"; exit 1 }
Ok "using .venv (nothing outside it is modified)"

Say "Installing dependencies"
& $VPy -m pip install --upgrade pip 2>&1 | Out-Null
& $VPy -m pip install -r requirements.txt 2>&1 | Out-File -FilePath "$env:TEMP\mosshoos_pip.log"
if ($LASTEXITCODE -eq 0) {
    Ok "numpy pandas scipy scikit-learn matplotlib"
} else {
    Bad "pip install failed - see $env:TEMP\mosshoos_pip.log"
    Get-Content "$env:TEMP\mosshoos_pip.log" -Tail 5
    exit 1
}

# --- 4. verify ---------------------------------------------------------------
Say "Verifying the method (unit tests)"
$Fail = $false
foreach ($t in @("sddp_clean", "sde_convergence", "hjb_clean", "wavelet_galerkin")) {
    Write-Host ("    {0,-18} " -f $t) -NoNewline
    $out = & $VPy "src\$t.py" 2>&1 | Out-String
    if ($out -match "ALL TESTS PASSED") {
        Write-Host "PASS" -ForegroundColor Green
    } else {
        Write-Host "FAIL" -ForegroundColor Red
        $Fail = $true
        ($out -split "`n" | Select-String "TEST .*(PASS|FAIL)") | ForEach-Object { Write-Host "           $_" }
    }
}
if (-not $Fail) { Ok "all components validated" } else { Bad "some tests failed (see above)" }

# --- 5. run ------------------------------------------------------------------
if (-not $NoRun) {
    Say "Running the pipeline on real OMVS data (about a minute)"
    & $VPy compare_plan4res.py --inflow data/OMVS_inflows.csv --n_iter 150 --n_mc 300 `
        2>&1 | Out-File -FilePath "$env:TEMP\mosshoos_run.log"
    if ($LASTEXITCODE -eq 0) { Ok "pipeline complete" }
    else { Bad "pipeline failed"; Get-Content "$env:TEMP\mosshoos_run.log" -Tail 8 }
    Get-Content "$env:TEMP\mosshoos_run.log" | Select-String "^  (MS|MS0|MD|MD0) " |
        ForEach-Object { Write-Host "    $_" }

    if (-not (Test-Path "results\storage_sweep.json")) {
        Say "Storage-sensitivity sweep (a few minutes)"
        & $VPy storage_sweep.py --inflow data/OMVS_inflows.csv --out results/storage_sweep.json `
            2>&1 | Out-File -FilePath "$env:TEMP\mosshoos_sweep.log"
        if ($LASTEXITCODE -eq 0) { Ok "sweep complete" } else { Bad "sweep failed" }
    }

    Say "Building the dashboard"
    & $VPy make_static.py 2>&1 | ForEach-Object { Write-Host "    $_" }
    $D = Join-Path (Get-Location) "dashboard.html"
    if (Test-Path $D) { Start-Process $D; Ok "opened $D" }
}

# --- done --------------------------------------------------------------------
$loc = (Get-Location).Path
Write-Host @"

  ------------------------------------------------------------------
  Done.  Everything below assumes you are in
  $loc
  and have run:   .\.venv\Scripts\Activate.ps1

    start dashboard.html                      standalone, offline, shareable
    python dashboard_server.py --port 8050    live: change params, re-run
                                              then open http://localhost:8050

    python compare_plan4res.py --inflow data/OMVS_inflows.csv --n_iter 300
    python storage_sweep.py    --inflow data/OMVS_inflows.csv
    python src\sddp_clean.py                  re-verify any component

  Results: results\comparison_results.json  -  RESULTS.md
  ------------------------------------------------------------------

"@
