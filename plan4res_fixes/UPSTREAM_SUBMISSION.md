# Proposed upstream fix: MILPSolver set_par overload hiding disables LP relaxation
Affected: SMS++ MILPSolver; breaks Plan4RES SSV (hence SIM/CEM) on HiGHS-only builds.
Symptom: SDDPOptimizer::oneStepBackward: no linearization is available -> Aborted.
Cause: set_par(idx_type,int) in MILPSolver hides base set_par(idx_type,double);
relax_int_vars is narrowed into the int index space, UC solved as MIP, no duals, no cuts.
Fix (one line): using CDASolver::set_par;  (patch in this directory).
Verified: HiGHS solves LP, SSV writes cuts ('Optimal solution found.'), SIM+CEM run.
Repro: https://github.com/belvinos91/mosshoos-plan4res (tag v1.0).
