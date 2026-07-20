# Modifications to stock Plan4RES / SMS++ (all in plan4res_fixes/)
1. MILPSolver.h: using CDASolver::set_par; - restores LP relaxation (core fix).
2. boost::mpi::communicator moved after MPI_Init (file-scope crash).
3. Added missing -DUSE_MPI compile definition.
4. ssv.sh: six bare-filename solver invocations -> full paths.
Effect: SSV->SIM->CEM run on HiGHS, no commercial licence. All else unmodified upstream.
