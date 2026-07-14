# Plan4RES / SMS++ sddp_solver — fixes applied

To reproduce the Plan4RES baseline we had to fix several issues in the
SMS++ sddp_solver build and the p4r wrapper scripts.

## 1. MPI initialization (sddp_solver.cpp)
File-scope `boost::mpi::communicator world` called MPI_Comm_rank before
MPI_Init. Moved communicator construction inside main() after
`boost::mpi::environment env(argc,argv)`; kept a file-scope pointer +
`#define world (*plan4res_world_ptr)` so the deferred-MPI helper still compiles.

## 2. USE_MPI compile flag (tools/sddp_solver/CMakeLists.txt)
The `#ifdef USE_MPI` block in main() was compiled out because -DUSE_MPI was
missing from CXX_DEFINES. Added:
    target_compile_definitions(sddp_solver PRIVATE USE_MPI)

## 3. p4r path passing (scripts/include/ssv.sh)
p4r invoked sddp_solver with bare filenames that only resolved from a specific
working directory. Replaced with full paths on all 6 solver-call lines:
  block : ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4
  config: -S ${CONFIG_IN_P4R}/sddp_solver.txt

## STATUS
Fixes 1–3 resolved MPI init, block loading, config loading. HiGHS solves the
UC subproblems successfully.

## REMAINING (reported to maintainers)
`SDDPOptimizer::oneStepBackward: no linearization is available` on the stock
toyDataset — the SDDPBlock has no last-stage cut and this SMS++ build aborts
rather than using the all-zero fallback. HiGHS is the only compiled MILP solver
(SCIP/Gurobi not built/licensed).

## WORKAROUND (this repository)
`src/plan4res_io.py` reads Plan4RES's own toyDataset inputs, and the MOSSHOOS
pipeline reproduces the Plan4RES deterministic method (M-D) directly — enabling
the M-S vs M-D comparison without depending on the unfixed solver.
