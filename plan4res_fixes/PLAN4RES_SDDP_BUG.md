# FIXED (patch below): SDDP fails on HiGHS-only builds - `set_par` overload hiding silently overwrites `intRelaxIntVars`

**Component:** SDDPBlock / SDDPSolver, MILPSolver (HiGHS)
**Severity:** blocks the **entire Plan4RES workflow** on open-source-solver builds
(SSV fails; SIM and CEM then fail because they consume SSV's Bellman values)
**Reproducible:** yes, on the stock `toyDataset` shipped with Plan4RES
**Status:** ROOT CAUSE FOUND, one-line fix verified - see `0001-fix-MILPSolver-overload-hiding.patch`

---

## Root cause (found, fixed, verified)

A C++ **overload-hiding** bug in `MILPSolver`.

`MILPSolver.h` declares:

```cpp
534:  void set_par( idx_type par , int value ) override;
537:  // void set_par( idx_type par , double value ) override;   <-- commented out
540:  void set_par( idx_type par , std::string && value ) override;
```

Declaring *any* `set_par` in a derived class **hides every base-class `set_par`
overload**. The `double` version is commented out and there is **no
using-declaration**, so `MILPSolver` has no `double` overload at all - neither
declared nor inherited.

`HiGHSMILPSolver::set_par( par , double )` ends with:

```cpp
MILPSolver::set_par( par , value );      // value is a double
```

With no `double` overload to bind to, the compiler **narrows the double to int**
and dispatches to `set_par( idx_type , int )`.

**The int and double parameter enums are separate index spaces.** So in the stock
`uc_solverconfig.txt`:

```
dblAAccSol 5e-04        # double parameter, index 7
```

`5e-04` truncates to `int 0`, calls `set_par( 7 , 0 )`, and int index 7 is
**`intRelaxIntVars`** - which is therefore silently set to false.

```
dblAAccSol 5e-04
  -> set_par( 7 , 0 )              (double narrowed to int)
  -> relax_int_vars = false        (int idx 7 == intRelaxIntVars)
  -> UC solved as a MIP            (204 integer variables, 187 binary)
  -> a MIP has no dual solution
  -> HiGHSMILPSolver::has_dual_solution() == false   (kHighsSolutionStatusNone)
  -> BendersBFunction::has_linearization( true ) == false
  -> SDDPSolver::SDDPOptimizer::oneStepBackward throws
```

This is invisible to inspection - every layer of the configuration path is
correct. It was found by instrumenting `MILPSolver::set_par` and taking a
backtrace on the offending call:

```
[DBG] set_par intRelaxIntVars <- 1  (this=0x7bd7380)      <- config applies it
[DBG] set_par(INT) par=7 value=0 intRelaxIntVars_idx=7    <- dblAAccSol lands here
[DBG] load_problem: relax_int_vars=0                      <- MIP

backtrace:
  HiGHSMILPSolver::set_par(unsigned short, double)         <- the double overload
  ThinComputeInterface::set_ComputeConfig
  BlockSolverConfig::apply(Block*)
  SDDPSolver::set_Block
  Block::register_Solver
  BlockSolverConfig::apply(Block*)
  process_block_file
  main
```

## The fix

One line in `MILPSolver/include/MILPSolver.h`:

```cpp
using CDASolver::set_par;          // restore the hidden base-class overloads
void set_par( idx_type par , int value ) override;
```

## Verified

Stock toyDataset, HiGHS 1.15.1, **no commercial solver**:

| | before | after |
|---|---|---|
| `relax_int_vars` at load | 0 | **1** |
| HiGHS problem type | `Solving MIP ... 204 integer variables` | **`Solving LP`** (203 LP solves, 0 MIP) |
| linearization errors | fatal throw | **0** |
| exit code | Aborted (core dumped) | **0, "Optimal solution found."** |
| `p4r SSV toyDataset` | fails | **completes, writes `results_simul/cuts.txt`** |

Backward value: 32818090192.099266052
Forward value:  37582251925.978507996

## Summary

On a build where **HiGHS is the only compiled MILP solver**, `sddp_solver` aborts on the
first backward step with:

```
terminate called after throwing an instance of 'std::logic_error'
  what():  SDDPOptimizer::oneStepBackward: no linearization is available.
```

The root cause is that the unit-commitment subproblem is solved as a **MIP** despite
`intRelaxIntVars 1` being set in `uc_solverconfig.txt`. A MIP has no dual solution, so no
linearization of the `BendersBFunction` can be produced, and SDDP cannot generate a cut.

This is not dataset-specific, and it is not limited to SSV. Because `SIM` and `CEM`
consume the Bellman values that SSV produces, the failure cascades:

```
p4r FORMAT toyDataset -M simul   ->  successfully ran FORMAT
p4r SIM    toyDataset            ->  "No Bellman values found in .../results_simul/
                                      nor .../results_optim/"
                                     "error while running SIM"
p4r FORMAT toyDataset -M invest  ->  successfully ran FORMAT
p4r CEM    toyDataset            ->  "None of BellmanValuesOUT.csv and cuts.txt is
                                      present in results_optim/"
                                     "error while running CEM with investment_solver"
```

So on a HiGHS-only build only `CREATE` and `FORMAT` succeed; **SSV, SIM, CEM and
SSVandCEM are all unusable**. It appears to affect every user without a commercial solver
licence, which we suspect is why it has not surfaced upstream: CPLEX/Gurobi/SCIP builds
presumably take a path where the relaxation is applied.

---

## Reproduce

```bash
cd ~/plan4res/p4r-env && ./bin/p4r
CFG=<instance>/settings
sed -i 's/^intLogVerb 0/intLogVerb 2/' $CFG/uc_solverconfig.txt   # show HiGHS's problem type

sddp_solver -S "$CFG/sddp_solver.txt" -c "$CFG/" \
  -p <instance>/nc4_optim/ <instance>/nc4_optim/SDDPBlock.nc4
```

## Observed

```
Warning: SDDPSolver::compute: No cut for the last stage has been provided and
         ... an all-zero cut will be used.

Solving MIP model with:
MIP has 757 rows; 541 cols; 1802 nonzeros; 204 integer variables (187 binary)
Presolve reductions: rows 276(-481); columns 310(-231); nonzeros 983(-819)
  Solution status   feasible
                    3009586843.29 (objective)

terminate called after throwing an instance of 'std::logic_error'
  what():  SDDPOptimizer::oneStepBackward: no linearization is available.
Aborted (core dumped)
```

**`Solving MIP model with: ... 204 integer variables (187 binary)`** is the problem.
With `intRelaxIntVars 1` the subproblem should be solved as an LP.

## Expected

The inner UCBlock is solved as an LP; HiGHS returns a dual solution; the
`BendersBFunction` produces a diagonal linearization; SDDP builds the cut.

---

## Mechanism (traced)

```
intRelaxIntVars 1 does not take effect on the inner UCBlock at solve time
  -> HiGHS solves a MIP                                  (204 integer variables)
  -> a MIP has no dual solution
  -> HiGHSMILPSolver::has_dual_solution()                (HiGHSMILPSolver.cpp:693)
       Highs_getIntInfoValue(highs,"dual_solution_status",&s)
       s == kHighsSolutionStatusNone  ->  return false
  -> BendersBFunction::has_linearization( true )  ->  false
  -> SDDPSolver::SDDPOptimizer::oneStepBackward           (SDDPSolver.cpp:1265)
       else { throw std::logic_error(
              "SDDPOptimizer::oneStepBackward: no linearization is available." ) }
```

## What we verified (all correct - the defect is elsewhere)

We checked every layer of the configuration path and found no fault in any of them:

| Layer | File / line | Status |
|---|---|---|
| Config file sets the parameter | `uc_solverconfig.txt` | `intRelaxIntVars 1` present |
| `strInnerBSC` is read | `SDDPSolver.h:2146` | OK |
| `-c` prefix is applied | `sddp_solver.cpp:263` `Configuration::set_filename_prefix` | OK |
| `ssv.sh` passes a trailing slash | `ssv.sh:114,120` `-c ${CONFIG_IN_P4R}/` | OK |
| Inner BSC is deserialized, errors are thrown not swallowed | `SDDPSolver.cpp:189-197` | OK |
| Inner BSC applied to **every** inner block | `SDDPSolver.cpp:~218` `f_inner_block_solver_config->apply( inner_block )` | OK |
| `apply()` configures **before** registering | `BlockSolverConfig.cpp` (both add and replace paths) | OK |
| String -> enum mapping | `MILPSolver.cpp:1957` `"intRelaxIntVars" -> intRelaxIntVars` | OK |
| Enum -> field | `MILPSolver.cpp:1838` `relax_int_vars = bool(value)` | OK |
| HiGHS override delegates to parent | `HiGHSMILPSolver.cpp:1872` `MILPSolver::set_par(par,value)` | OK |
| HiGHS honours the flag when loading vars | `HiGHSMILPSolver.cpp:1032,1709` `if( var->is_integer() && !relax_int_vars )` | OK |
| Reset-to-defaults would need `f_diff == 0` | both ComputeConfigs declare `1` | not triggered |

Every layer reads correctly, yet the integers survive into the solve. We could not
identify where the flag is lost without instrumenting a rebuild, and we expect the
maintainers will recognise the initialization-order issue immediately.

## Both solver paths fail identically

`SDDPGreedySolver` configures the inner blocks by calling `configure_inner_block(stage)`
**explicitly at solve time** (`SDDPGreedySolver.cpp:89,214,326`), whereas `SDDPSolver`
relies on the `set_Block` path. We tested whether the two differ:

```bash
sddp_solver -s -S settings/sddp_greedy.txt -c settings/ -p nc4_optim/ nc4_optim/SDDPBlock.nc4
```

```
MIP has 756 rows; 541 cols; 1801 nonzeros; 204 integer variables (187 binary)
The subproblem at stage 0 is unbounded.
Lower bound: -inf
Upper bound: inf
```

**The greedy solver sees the same 204 integer variables.** The relaxation therefore fails
through *both* configuration paths, which rules out a path-specific plumbing bug and
points at something common to both - `BlockSolverConfig::apply()` itself, or the point at
which `MILPSolver` / `HiGHSMILPSolver` loads variables relative to when `relax_int_vars`
is set.

Note also the greedy run reports **"The subproblem at stage 0 is unbounded"**: with an
all-zero last-stage cut and no relaxation, the stage problem has no finite bound.

## Suggested diagnostics

1. Log `relax_int_vars` inside `HiGHSMILPSolver` immediately before the variables are
   loaded (`HiGHSMILPSolver.cpp:1032,1709`). We expect it to be `false` there despite
   `intRelaxIntVars 1`.
2. Confirm whether the inner `Solver` instance configured by
   `f_inner_block_solver_config->apply( inner_block )` is the same instance later
   returned by `BendersBFunction::get_solver()`.
3. Check whether `set_par( intRelaxIntVars , 1 )` arrives **after** the block's variables
   have already been loaded into HiGHS. `MILPSolver::set_par` only sets a flag
   (`MILPSolver.cpp:1838-1839`); nothing re-reads it, so a model already built with
   integers would keep them.

## Secondary issues noticed

- **Silent skip.** `sddp_solver.cpp:1260-1264` does
  `std::ifstream f( config_filename_prefix + strInnerBSC ); if( !f.is_open() ) continue;`
  - a missing inner config is silently ignored. A warning would make failures like this
  far easier to diagnose. Note also the bare `+` concatenation: a `-c` argument without a
  trailing separator silently produces an invalid path.
- **Last-stage cut.** `SDDPSolver::compute` warns "No cut for the last stage has been
  provided" on the stock toyDataset, then proceeds with an all-zero cut. Whether the
  shipped dataset should include a terminal cut may be worth confirming separately.

## Environment

- HiGHS 1.15.1 (git 04024d701f)
- MILP solvers compiled: **HiGHS only** (no CPLEX / Gurobi / SCIP)
- Dataset: stock `toyDataset`, RoGonDor reservoir, `p4r FORMAT ... -M optim`
- Invocation: `p4r SSV`, `p4r SSVandCEM`, and direct `sddp_solver` - identical failure
- Unaffected by `NumberSSVIterations` (3 -> 20) or `NumberSSVForwardFirstStep` (1 -> 4)

## Related fixes we needed to get this far

Three further defects blocked `sddp_solver` before reaching this point; patches available
if useful:

1. **MPI initialization.** A file-scope `boost::mpi::communicator world` calls
   `MPI_Comm_rank` before `MPI_Init` (which happens in `main` via
   `boost::mpi::environment`). Moving the communicator into `main` after `env(argc,argv)`
   fixes the abort.
2. **`USE_MPI` never defined.** The `#ifdef USE_MPI` block in `main` was compiled out;
   `target_compile_definitions(sddp_solver PRIVATE USE_MPI)` was missing from
   `tools/sddp_solver/CMakeLists.txt`.
3. **Bare filenames in `ssv.sh`.** Six solver-invocation lines passed `SDDPBlock.nc4` and
   `sddp_solver.txt` without paths, resolving only from one specific working directory.
