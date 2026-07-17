/*--------------------------------------------------------------------------*/
/*---------------------- File SDDPGreedySolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the SDDPGreedySolver class, implementing the Solver
 * interface, for multistage programming problems defined by the
 * SDDPBlock. The SDDPGreedySolver implements a greedy strategy to try to
 * solve a deterministic (single-scenario) multistage problem encoded by an
 * SDDPBlock as defined below.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright Copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SDDPGreedySolver
#define __SDDPGreedySolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SDDPBlock.h"
#include "Solver.h"

#include <random>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class BendersBFunction;      // forward declaration of BendersBFunction

 class BlockSolverConfig;     // forward declaration of BlockSolverConfig

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SDDPGreedySolver_CLASSES Classes in SDDPGreedySolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS SDDPGreedySolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a greedy solver for multistage programming problems
/**
 * The SDDPGreedySolver class derives from Solver and implements a sequential,
 * greedy strategy to solve an SDDPBlock for a fixed scenario. Recall that an
 * SDDPBlock represents an optimization problem of the form
 *
 * \f[
 *   \min_{x_0 \in \mathcal{X}_0} f_0(x_0) +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_1 \in \mathcal{X}_1} f_1(x_1) +
 *   \mathbb{E} \left \lbrack \dots +
 *   \mathbb{E} \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{X}_{T-1}} f_{T-1}(x_{T-1})
 *   \right\rbrack \right\rbrack\right\rbrack, \qquad (1)
 * \f]
 *
 * where T is the time horizon, \f$\mathcal{X}_t \equiv
 * \mathcal{X}_t(x_{t-1}, \xi_t) \subseteq \mathbb{R}^{n_t}\f$ for each
 * \f$t \in \{0, \dots, T-1\}\f$, and \f$ \xi = \{ \xi_t \}_{t \in \{1, \dots,
 * T-1\}} \f$ is a stochastic process. See SDDPBlock for details. The
 * SDDPGreedySolver considers the problem (1) for a single realization of the
 * stochastic process, i.e., a deterministic problem of the form
 *
 * \f[
 *   \min_{x_0 \in \mathcal{\tilde{X}}_0} f_0(x_0) +
 *   \left \lbrack
 *   \min_{x_1 \in \mathcal{\tilde{X}}_1} f_1(x_1) +
 *   \left \lbrack \dots +
 *   \left \lbrack
 *   \min_{x_{T-1} \in \mathcal{\tilde{X}}_{T-1}} f_{T-1}(x_{T-1})
 *   \right\rbrack \right\rbrack\right\rbrack, \qquad (2)
 * \f]
 *
 * with \f$\mathcal{\tilde{X}}_t \equiv \mathcal{\tilde{X}}_t(x_{t-1},
 * \tilde{\xi}_t)\f$ where \f$ \tilde{\xi}_t = \{ \tilde{\xi}_t \}_{t \in \{1,
 * \dots, T-1\}} \f$ is a realization of the stochastic process \f$ \xi
 * \f$. The SDDPGreedySolver is a heuristic as it does not look for an optimal
 * solution to problem (2). The method it employs is very simple: it solves
 * the subproblem at each stage in sequence, from the first to the last one,
 * using the solution found for a stage to define the problem at the next
 * stage. First, for a given \f$ (x_{-1}, \tilde{\xi}_0)\f$, it solves
 * the problem
 *
 * @f{align}
 *   \min       & \ \ f_0(x_0) + \mathcal{P}_{1}(x_0) \qquad (3) \\
 *   {\rm s.t.} & \ \ x_0 \in \mathcal{\tilde{X}}_0(x_{-1},
 *                            \tilde{\xi}_0)
 * @f}
 *
 * where \f$ \mathcal{P}_{1} \f$ is a polyhedral function that approximates
 * the cost-to-go function. Let \f$ x_0^* \f$ be a solution obtained to
 * problem (3). Next, the following problem is solved
 *
 * @f{align}
 *   \min       & \ \ f_1(x_1) + \mathcal{P}_{2}(x_1)\\
 *   {\rm s.t.} & \ \ x_1 \in \mathcal{\tilde{X}}_1(x^*_{0},
 *                  \tilde{\xi}_1)
 * @f}
 *
 * and a solution \f$ x_1^* \f$ is obtained. This process continues until the
 * subproblem at stage \f$ T-1 \f$ is solved and a solution \f$ x_{T-1}^* \f$
 * is found for it. In general, for each \f$ t \in \{0, \dots, T-1\}\f$, the
 * subproblem solved at stage \f$ t \f$ is the following
 *
 * @f{align}
 *   \min       & \ \ f_t(x_t) + \mathcal{P}_{t+1}(x_t)\\
 *   {\rm s.t.} & \ \ x_t \in \mathcal{\tilde{X}}_t(x^*_{t-1},
 *                  \tilde{\xi}_t)
 * @f}
 *
 * where \f$ x^*_{t-1} \f$ is a solution to the subproblem at stage \f$ t-1
 * \f$ if \f$ t > 0 \f$, \f$ x^*_{-1} \equiv x_{-1} \f$, and \f$
 * \mathcal{P}_{t+1} \f$ denotes a polyhedral function that approximates the
 * cost-to-go function at stage \f$ t \f$.
 *
 *     Notice that the SDDPGreedySolver does not solve neither the problem
 *     encoded by SDDPBlock nor the deterministic (single-scenario) multistage
 *     problem defined in (2). It may not even find a feasible solution to
 *     problem (2) even if one exists.
 */

class SDDPGreedySolver : public CDASolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 /// "import" Index from SDDPBlock
 using Index = SDDPBlock::Index;

 /// public enum for the possible return values of compute()
 /** Public enum "extending" Solver::sol_type with more detailed values
  * specific to SDDPGreedySolver. Among the values defined in
  * Solver::sol_type, only kOK is not considered in SDDPGreedySolver. A few
  * comments about the specific interpretation of these values is in order. As
  * explained above, the SDDPGreedySolver solves each subproblem
  * sequentially. The following values are returned by compute() depending on
  * the value returned by compute() when solving each subproblem.
  *
  * - kError is returned when an unrecoverable error occurs while solving a
  *   subproblem. The stage of the subproblem at which the error occurred can
  *   then be retrieved by the method get_fault_stage().
  *
  * - kUnbounded is returned when some subproblem is unbounded. The stage of
  *   the unbounded subproblem can then be retrieved by the method
  *   get_fault_stage(). It does not mean, however, that the subproblems at
  *   later stages are all feasible.
  *
  * - kInfeasible is returned when the subproblem at the first stage is
  *   infeasible. This is the only case in which the SDDPGreedySolver
  *   guarantees that the deterministic (single-scenario) multistage problem
  *   is infeasible. In this case, instead of returning kInfeasible we could
  *   return kSubproblemInfeasible (see below), but we decided to return
  *   kInfeasible as this value should be more broadly understood.
  *
  * - kStopTime is returned when all subproblems were "solved", but the
  *   solution of some subproblem terminated with a kStopTime status and all
  *   subproblems at previous stages terminated with either a kOK or a
  *   kLowPrecision status. The stage at which this event occurred can then be
  *   retrieved by the method get_fault_stage(). In other words,
  *   get_fault_stage() will return t such that the solution of subproblem at
  *   stage t terminated with a kStopTime status and the solution of each
  *   subproblem at stage in {0, ..., t-1} terminated with either a kOK or a
  *   kLowPrecision status.
  *
  * - kStopIter is returned when all subproblems were "solved", but the
  *   solution of some subproblem terminated with a kStopIter status and all
  *   subproblems at previous stages terminated with either a kOK or a
  *   kLowPrecision status. The stage at which this event occurred can then be
  *   retrieved by the method get_fault_stage(). In other words,
  *   get_fault_stage() will return t such that the solution of subproblem at
  *   stage t terminated with a kStopIter status and the solution of each
  *   subproblem at stage in {0, ..., t-1} terminated with either a kOK or a
  *   kLowPrecision status.
  *
  * - kLowPrecision is returned when every subproblem is "solved" and
  *   terminated with either a kOK or kLowPrecision status. In this case, the
  *   solution found is feasible for the deterministic (single-scenario)
  *   multistage problem but there is no guarantee that it is optimal.
  */

 enum sddp_greedy_sol_type {
 kSubproblemInfeasible = kInfeasible + 1 ,
 ///< some subproblem may be infeasible
 /**< It means that a subproblem at some stage, let say t, different than the
  * first one turned out to be infeasible. Since the feasible region of a
  * subproblem may depend on the solution of the subproblem at the previous
  * stage, it does not mean that the deterministic (single-scenario)
  * multistage problem is infeasible, as there could be another solution for
  * the subproblem at the previous stage that would make the subproblem at
  * stage t feasible. The stage t of the infeasible problem can be retrieved
  * by the method get_fault_stage().
  */

 kSolutionNotFound ,
 ///< the solution to some subproblem has not been found
 /**< It means that a solution to a subproblem has not been found for whatever
  * reason. The stage at which the solution could not be found can be
  * retrieved by the method get_fault_stage().
  */

 };  // end( sddp_greedy_sol_type )

/*--------------------------------------------------------------------------*/

 /// public enum for the int algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters of
  * "int" type that the SDDPGreedySolver has, besides those defined in
  * Solver. The value intLastAlgPar is provided so that the list can be easily
  * further extended by derived classes. */

 enum int_par_type_SDDP_Greedy_S {

  intScenarioId = int_par_type_S::intLastAlgPar ,
  ///< The id of the scenario that must be considered
  /**< This is the id of the scenario that must be considered when trying to
   * solve the deterministic (single-scenario) multistage problem. It must be
   * a valid id for a scenario handled by the SDDPBlock. By default, its value
   * is 0, which means that the scenario with id 0 will be considered when
   * solving the subproblems at each stage (except possibly at the first
   * stage; see #intFirstStageScenarioId). */

  intFirstStageScenarioId ,
  ///< The id of the scenario to be considered at the first stage
  /**< This parameter specifies the id of the scenario that must be considered
   * while solving the subproblem at the first stage. If it is -1, then the id
   * of the scenario to be considered at the first stage is that specified by
   * #intScenarioId. Otherwise, if it is negative (less than -1), it means
   * that no scenario must be set while solving the sub-problem at the first
   * stage (i.e., the data for that subproblem has already been set, except
   * possibly the initial state). If it is nonnegative, then it is the id of
   * the scenario to be considered at the first stage and thus it must be a
   * number between 0 and the total number of scenarios minus 1. By default,
   * its value is -1 (i.e., the id of the scenario for the first stage is that
   * given by #intScenarioId). */

  intUnregisterSolver ,
  ///< Indicates whether Solver(s) of the inner Block must be unregistered
  /**< Within compute(), for each stage, the inner Block of the
   * BendersBFunction at that stage is solved. Sometimes, the Solver(s)
   * attached to that inner Block are not needed after it is solved and can be
   * unregistered. This parameter indicates whether these Solver(s) must be
   * unregistered. If the value of this parameter is non-zero, then the
   * Solver(s) of the inner Block of the BendersBFunction are unregistered as
   * follows. If this SDDPGreedySolver has a BlockSolverConfig associated with
   * that inner Block, then it is used to unregister the Solver(s) of that
   * inner Block. Otherwise, all Solver(s) of that inner Block are
   * unregistered by a call to Block::unregister_Solvers(). */

  intScenarioSeed ,
  ///< Seed for the random number engine that selects the scenarios
  /**< This parameter determines whether the scenarios to be considered at
   * each stage (from the second stage onwards) should be randomly selected.
   * The id of the scenario for the first stage problem is always determined
   * by #intFirstStageScenarioId. If #intScenarioSeed is negative, then the
   * scenarios to be considered are determined by #intScenarioId (and
   * #intFirstStageScenarioId). If #intScenarioSeed is nonnegative, then the
   * scenarios for all stages except the first one are selected at random and
   * #intScenarioSeed serves as the seed for the random number engine that
   * selects the scenarios. The default value for this parameter is negative,
   * which means that the scenarios are not selected at random and, thus, are
   * determined by #intScenarioId (and #intFirstStageScenarioId). Notice that
   * this parameter has a higher priority over the parameter
   * #intScenarioId. This means that if #intScenarioSeed is nonnegative, then
   * the scenarios for each stage from the second stage onwards are selected
   * at random, no matter the value of #intScenarioId. */

  intScenarioSampleFrequency ,
  ///< Frequency at which scenarios should be sampled
  /**< This parameter determines the frequency at which scenarios should be
   * sampled when the scenarios are chosen at random (see the #intScenarioSeed
   * parameter). If it is positive, scenarios are sampled every
   * #intScenarioSampleFrequency stages. When a scenario is not sampled for
   * some stage, the id of the scenario for that stage will be the same as the
   * id of the scenario for the previous stage. If the value for this
   * parameter is nonpositive, then scenarios are not sampled at all and,
   * therefore, the id of the scenario for each stage will be the same as the
   * id of the scenario for the first stage. The default value for this
   * parameter is 1, which means that scenarios are sampled for every
   * stage. */

  intSimulationDataOutputPrecision ,
  ///< Precision of the simulation data output
  /**< This parameter determines the precision of the simulation data that is
   * output (see #strSimulationData). The value of this parameter will be
   * passed to the std::setprecision() function when the simulation data is
   * output. The default value for this parameter is 20. */

  intOutputScenario ,
  ///< It indicates how (and if) scenarios must be output
  /**< This parameter determines how (and if) the scenario considered during
   * the simulation must be output if the parameter #strSimulationData is
   * provided. If #intOutputScenario is negative, then only the IDs of the
   * scenarios considered at each stage are output. If it is positive, then
   * the scenarios themselves are output. If it is zero, then no scenario
   * (neither the ID nor the scenario itself) is output. By default, the value
   * of this parameter is -1, which means that the IDs of the scenarios are
   * output (if #strSimulationData is provided). */

  intEarlyConfig ,
  ///< It indicates whether the inner Blocks should be configured in advance
  /**< By default, if a BlockConfig or a BlockSolverConfig has been provided
   * (by means of the parameters #strInnerBC and #strInnerBSC, respectively,
   * or by the "extra" Configuration in the ComputeConfig of this
   * SDDPGreedySolver, see set_ComputeConfig()), the inner Block of each
   * BendersBFunction is configured right before it is solved, within
   * SDDPGreedySolver::compute(). However, in some situations it may be
   * necessary to configure every inner Block in advance, at the time this
   * SDDPGreedySolver is configured. This parameter determines when the inner
   * Blocks are configured. If #intEarlyConfig is nonzero, then the inner
   * Blocks are configured sooner, at the time the SDDPGreedySolver is
   * configured. By default, #intEarlyConfig is zero, which means that each
   * inner Block is configured right before it is solved. */

  intLoadCutsOnce ,
  ///< It indicates whether cuts should be loaded only once for each stage
  /**< This parameter specifies how many times the cuts provided by
   * #intLoadCuts are loaded. Cuts associated with a particular stage, if
   * provided, are loaded within compute() right before the subproblem
   * associated with that stage is solved. If this parameter is nonzero, the
   * cuts are loaded only once. That is, if cuts associated with a particular
   * stage are loaded within a call to compute(), then no cuts are loaded for
   * that stage in any subsequent calls to compute(). If this parameter is
   * zero, then cuts are loaded at each call to compute().
   *
   * Notice that cuts are not necessarily loaded during the first call to
   * compute(). In a normal situation, when all subproblems are solved within
   * a call to compute(), then all cuts are loaded in that call. However, not
   * all subproblems may be solved within a call to compute(), which may
   * happen, for instance, when the subproblem associated with a certain stage
   * S is infeasible. In this case, no cuts are loaded for any stage T greater
   * than S and any cuts associated with T may only be loaded in a possible
   * next call to compute().
   *
   * By default, #intLoadCutsOnce is 1, which means that cuts are loaded only
   * once. */

  intLastAlgPar
  ///< first allowed new double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set of
   * int algorithmic parameters. */

 };  // end( int_par_type_SDDP_Greedy_S )

/*--------------------------------------------------------------------------*/

 /// public enum for the string algorithmic parameters
 /** Public enum describing the different types of algorithmic parameters of
  * "string" type that the SDDPGreedySolver has, besides those defined in
  * Solver. The value strLastAlgPar is provided so that the list can be easily
  * further extended by derived classes. */

 enum str_par_type_SDDP_Greedy_S {

  strInnerBC = str_par_type_S::strLastAlgPar ,
  ///< name of the file containing the default BlockConfig for the inner Block
  /**< Name of the file containing the default BlockConfig that will be
   * applied to the inner Block of each BendersBFunction. */

  strInnerBSC ,
  ///< name of the file containing the default BlockSolverConfig for inner Block
  /**< Name of the file containing the default BlockSolverConfig that will be
   * applied to the inner Block of each BendersBFunction. */

  strLoadCuts ,
  ///< name of the file out of which cuts will be loaded
  /**< This parameter indicates the path to the file out of which cuts will be
   * loaded. By default, the path to this file is empty, which means that no
   * cut is loaded. If provided, the file must have the format the following
   * format. The first line contains the header, which will be simply
   * ignored. Each of the following lines must contain a cut described as
   * follows:
   *
   *     s, a_0, a_1, ..., a_{k-1}, b
   *
   * where s is a stage between 0 and get_time_horizon() - 1, which indicates
   * the stage with which the cut is associated, a_0, ..., a_{k-1} are the
   * coefficients of the cut (a_i being the coefficient associated with the
   * i-th state variable), and b is the constant (independent) term of the
   * cut. Cuts associated with a a particular stage are loaded within
   * compute() right before the subproblem associated with that stage is
   * solved. See the parameter #intLoadCutsOnce to control when cuts are
   * loaded. */

  strRandomCutsFile ,
  ///< name of the file out of which the random cuts will be retrieved
  /**< A random cut is a cut associated with a particular scenario. This
   * parameter indicates the path to the file out of which the random cuts
   * will be deserialized. By default, the path to this file is empty, which
   * means that no random cut is considered. If provided, the file must have
   * the format specified by SDDPBlock::deserialize_random_cuts(). */

  strSimulationData ,
  ///< name of the file in which data from the simulation will be saved
  /**< This is the name of the file in which some data obtained during the
   * simulation (i.e., the most recent call to compuet()), which includes
   * subgradients and values of the objectives of the subproblems, initial
   * states, and scenarios, will be saved. At each stage (except the first
   * one), the subgradients of the objective function with respect to the
   * initial state and with respect to the solution (final state) of the
   * subproblem at that stage are computed. At the first stage, only the
   * subgradient with respect to the solution of the first stage subproblem is
   * computed. At the end of a call to compute(), these subgradients are
   * output to the file whose name (path) is given by #strSimulationData. This
   * file will have the following format. The first line contains two
   * integers: the time horizon T and the number of initial states that will
   * be output (which is either T or T+1) separated by comma. This line is
   * followed by T or T+1 lines (depending on whether an initial state for the
   * first stage subproblem has been provided), each one containing the
   * initial state of some stage. Each of these lines have the following
   * format:
   *
   *     t, s_0, s_1, ..., s_{k-1}
   *
   * where t is a stage between 0 and T and (s_0, ..., s_{k-1}) is the initial
   * state (which has size k) for stage t. If an initial state has been
   * provided (see #vdblInitialState) or the SDDPBlock contains an initial
   * state (as returned by SDDPBlock::get_initial_state()), then T+1 initial
   * states are output (for each t in {0, ..., T}). Otherwise, T initial
   * states are output (for each t in {1, ..., T}). Notice that, although the
   * stages that we consider are 0, ..., T-1, an initial state for stage T
   * (which is a stage that has not been defined) is output. This is just the
   * final state (solution) of the last stage subproblem.
   *
   * After that, each of the next 2*T - 1 lines contains a subgradient of the
   * objective of the subproblem associated with a particular stage t and with
   * respect to either the initial state or the final state. There are T
   * subgradients with respect to the final state (one for each t in {0, ...,
   * T-1}) and T-1 subgradients with respect to the initial state (one for
   * each t in {1, ..., T-1}). Each of these lines has the following format:
   *
   *     t, s, a_0, ..., a_{k-1}
   *
   * where t is a stage between 0 and T-1, s is either the character I (for
   * initial) or the character F (for final) indicating that the subgradient
   * is with respect to the initial or the final state, respectively, and a_0,
   * ..., a_{k-1} are the elements of the subgradient.
   *
   * Then, each of the next T lines contains two elements
   *
   *     t, c
   *
   * where t is a stage between 0 and T-1 and c is the value of the objective
   * function of the subproblem associated with stage t, disregarding the
   * value of the cost-to-go function (or value function, future value
   * function, future cost function). That is, the objective value of the
   * subroblem is f = c + F, where F is the value of the cost-to-go function.
   *
   * Finally, depending on the value of the parameter #intOutputScenario, the
   * scenario considered during compute() (which can be set by the
   * #intScenarioId parameter or by the function set_scenario_id()) is output
   * in the last T lines. If #intOutputScenario is positive, then each of
   * these lines has the form
   *
   *     t, scenario_t
   *
   * where t is a stage in {0, ..., T-1} and scenario_t is a vector containing
   * the scenario for the stage t. If #intOutputScenario is negative, then
   * each of these lines has the form
   *
   *     t, scenario_t_id
   *
   * where t is a stage in {0, ..., T-1} and scenario_t_id is the ID of the
   * scenario considered for the stage t. If #intOutputScenario is zero, then
   * these lines are not output.
   *
   * By default, the path to this file is empty, which means that no data
   * obtained during the simulation will be output. */

  strLastAlgPar
  ///< first allowed new string parameter for derived classes
  /**< Convenience value for easily allow derived classes
   * to extend the set of string algorithmic parameters. */

 };  // end( str_par_type_SDDP_Greedy_S )

/*--------------------------------------------------------------------------*/

 /// public enum for the vector-of-int parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-int type that SDDPGreedySolver has in addition to these of
  * Solver. The value vintLastAlgPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum vint_par_type_SDDP_Greedy_S {

  vintStagesSample = vint_par_type_S::vintLastAlgPar ,
  ///< Stages at which a new scenario must be sampled
  /**< When scenarios are selected at random (see #intScenarioSeed), this
   * parameter specifies the stages at which scenarios must be sampled. If, at
   * a stage t > 0, a scenario is not sampled, then the id of the scenario to
   * be considered at stage t is the same as the id of the scenario that was
   * considered at stage t-1. Thus, if #vintStagesSample contains the stages S
   * = {t_0, t_1, ..., t_k}, with t_i < t_{i+1} for all i in {0, ..., k-1},
   * then scenarios are sampled at each stage t_i for i in {0, ..., k} and the
   * id to be considered at a stage t > 0 that does not belong to S is the
   * same as the id of the scenario that was (sampled and) considered at stage
   *
   *     t_j = min{ t_i in S | t >= t_i }
   *
   * Recall that the id of the scenario to be considered at the first stage
   * (stage 0) is determined by #intScenarioId. However, if 0 belongs to
   * #vintStagesSample, then the id of the scenario for the first stage is
   * also sampled and #intScenarioId is ignored.
   *
   * The #vintStagesSample parameter is an alternative to the parameter
   * #intScenarioSampleFrequency. The parameters #vintStagesSample and
   * #intScenarioSampleFrequency can be used together so that a scenario is
   * sampled at a stage t if t belongs to #vintStagesSample or t satisfies the
   * condition specified by #intScenarioSampleFrequency.
   *
   * Each element of this vector must be between 0 and get_time_horizon() -
   * 1. By default, this vector is empty. */

  vintLastAlgPar
  ///< first allowed new vector-of-int parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set of
   * vector-of-int parameters. */

 };  // end( vint_par_type_SDDP_Greedy_S )

/*--------------------------------------------------------------------------*/

 /// public enum for the vector-of-double parameters
 /** Public enum describing the different algorithmic parameters of
  * vector-of-double type that SDDPGreedySolver has in addition to these of
  * Solver. The value vdblLastAlgPar is provided so that the list can be
  * easily further extended by derived classes. */

 enum vdbl_par_type_SDDP_Greedy_S {

  vdblInitialState = vdbl_par_type_S::vdblLastAlgPar ,
  ///< the initial state for the first stage problem
  /**< The parameter for setting the initial state, i.e., the initial state to
   * be considered in the subproblem of the first stage. The size of this
   * vector must be equal to the size of the initial state and the i-th
   * element in this vector will be the value of the i-th initial state
   * variable of the first stage subproblem. If this vector is empty, no
   * initial state is set for the first stage problem. By default, this vector
   * is empty. */

  vdblLastAlgPar
  ///< first allowed new vector-of-double parameter for derived classes
  /**< Convenience value for easily allow derived classes to extend the set of
   * vector-of-double parameters. */

 };  // end( vdbl_par_type_SDDP_Greedy_S )

/**@} ----------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING SDDPGreedySolver --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing SDDPGreedySolver
 *  @{ */

 /// constructor
 SDDPGreedySolver( void ) { }

/*--------------------------------------------------------------------------*/

 /// destructor
 virtual ~SDDPGreedySolver();

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/

 /// set a given integer (int) numerical parameter
 /** Set a given integer (int) numerical parameter. Besides
  * considering the integer parameters defined in #int_par_type_S,
  * this function also accepts the following parameters:
  *
  * - #intScenarioId [0]
  *
  * - #intFirstStageScenarioId [-1]
  *
  * - #intUnregisterSolver [0]
  *
  * - #intScenarioSeed [-1]
  *
  * - #intScenarioSampleFrequency [1]
  *
  * - #intSimulationDataOutputPrecision [20]
  *
  * - #intOutputScenario [-1]
  *
  * - #intEarlyConfig [0]
  *
  * - #intLoadCutsOnce [1]
  *
  * Please refer to the #int_par_type_SDDP_Greedy_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 /* FIX (overload hiding): declaring any set_par() here hides ALL base-class
  * set_par() overloads. With no double overload and no using-declaration, a
  * set_par( par , <double> ) call narrows the double to int and dispatches to
  * the int overload - whose enum is a SEPARATE index space - silently
  * corrupting an unrelated parameter. Same defect as MILPSolver, where
  * dblAAccSol (double idx 7) overwrote intRelaxIntVars (int idx 7). */
 using CDASolver::set_par;
 void set_par( const idx_type par , const int value ) override {
  switch( par ) {
   case( intScenarioId ): set_scenario_id( value ); return;
   case( intFirstStageScenarioId ):
    set_first_stage_scenario_id( value ); return;
   case( intUnregisterSolver ): f_unregister_solver = value; return;
   case( intLogVerb ): log_verbosity = value; return;
   case( intScenarioSeed ): {
    if( value >= 0 ) {
     f_seed = value;
     random_number_engine.seed( f_seed );
    }
    else
     f_seed = Inf< Index >();
   }
   case( intScenarioSampleFrequency ):
    f_scenario_sample_frequency = value;
    return;
   case( intSimulationDataOutputPrecision ):
    f_simulation_data_output_precision = value;
    return;
   case( intOutputScenario ): f_output_scenario = value; return;
   case( intEarlyConfig ): f_early_config = value; return;
   case( intLoadCutsOnce ): f_load_cuts_once = value; return;
  }
  Solver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

 /// set a given string parameter
 /** Set a given string parameter. Set a given string parameter. Besides
  * considering the integer parameters defined in #str_par_type_S, this
  * function also accepts the following parameters:
  *
  * - #strInnerBC [""]: the filename of the "default" BlockConfig of the inner
  *   Block of the BendersBFunction(s). If non-empty(), this parameter is used
  *   to create a BlockConfig that is apply()-ed to the inner Block of all
  *   BendersBFunction. If left empty(), no BlockConfig is apply()-ed.
  *
  * - #strInnerBSC [""]: the filename of the "default" BlockSolverConfig of
  *   the inner Block of the BendersBFunction(s). If non-empty(), this
  *   parameter is used to create a BlockSolverConfig that is apply()-ed to
  *   the inner Block of all BendersBFunction. If left empty(), no
  *   BlockSolverConfig is apply()-ed. Note that each BendersBFunction does
  *   require a working Solver attached to its inner Block (unless the inner
  *   Solver can avoid it for some specially structured inner Block), so this
  *   will have to be provided in some way (for instance, it can come from the
  *   "extra" Configuration in the ComputeConfig of SDDPGreedySolver, see
  *   set_ComputeConfig()).
  *
  * - #strLoadCuts [""]: the filename of (path to) the file out of which cuts
  *   will be loaded. By default, the path to this file is empty, which means
  *   that no cut is loaded.
  *
  * - #strRandomCutsFile [""]: the filename of (path to) the file containing
  *   the random cuts (cuts associated with a particular scenario). By
  *   default, the path to this file is empty, which means that no random cut
  *   is considered. If provided, the file must have the format specified by
  *   SDDPBlock::deserialize_random_cuts().
  *
  * - #strSimulationData [""]: the filename of (path to) the file to which the
  *   data obtained during the simulation (if any) will be output. By default,
  *   the path to this file is empty, which means that no data is output.
  *
  * Please refer to the #str_par_type_SDDP_Greedy_S enumeration for a detailed
  * description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( const idx_type par , const std::string & value ) override {
  switch( par ) {
   case( strInnerBC ):
    f_inner_block_config_filename = value;
    return;
   case( strInnerBSC ):
    f_inner_block_solver_config_filename = value;
    return;
   case( strLoadCuts ):
    f_load_cuts_filename = value;
    return;
   case( strRandomCutsFile ):
    f_random_cuts_filename = value;
    return;
   case( strSimulationData ):
    f_simulation_data_filename = value;
    return;
  }
  Solver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

 /// set the vector-of-int paramaters of SDDPGreedySolver
 /** Set a given vector-of-int paramater. Besides considering the
  * vector-of-int parameters defined in #vint_par_type_S, this function
  * also accepts the following parameters:
  *
  * - #vintStagesSample
  *
  * Please refer to the #vint_par_type_SDDP_Greedy_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( idx_type par , std::vector< int > && value ) override {
  switch( par ) {
   case( vintStagesSample ): v_stages_to_sample = value; return;
  }
  Solver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

 /// set the vector-of-double paramaters of SDDPGreedySolver
 /** Set a given vector-of-double paramater. Besides considering the
  * vector-of-double parameters defined in #vdbl_par_type_S, this function
  * also accepts the following parameters:
  *
  * - #vdblInitialState
  *
  * Please refer to the #vdbl_par_type_SDDP_Greedy_S enumeration for a
  * detailed description of each of them.
  *
  * @param par A parameter to be set.
  *
  * @param value The value for the given parameter.
  */

 void set_par( idx_type par , std::vector< double > && value ) override {
  switch( par ) {
   case( vdblInitialState ): v_initial_state = value; return;
  }
  Solver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

 /// set the whole set of parameters of this SDDPGreedySolver in one blow
 /** This method sets the whole set of parameters of this SDDPGreedySolver in
  * one blow using a ComputeConfig object.
  *
  * Besides considering all the parameters of an SDDPGreedySolver, it can also
  * be used to configure the inner Block of every BendersBFunction by means of
  * the extra Configuration (ComputeConfig::f_extra_Configuration). If the
  * pointer to the extra Configuration is not nullptr, it can be any of the
  * following:
  *
  * - a pointer to a BlockConfig, which will be used to configure the inner
  *   Block of the BendersBFunction at every stage;
  *
  * - a pointer to a BlockSolverConfig, which will be used to configure the
  *   Solver of the inner Block of the BendersBFunction at every stage;
  *
  * - a pointer to a SimpleConfiguration< std::vector< Configuration * > >.
  *
  * In the last case, the first element of the vector, if present, must be
  * either nullptr or a pointer to a BlockConfig. The second element, if
  * present, must be either nullptr or a pointer to a BlockSolverConfig. These
  * will be used to configure the inner Block of the BendersBFunction at every
  * stage and their Solver. The third element, if present, must be either
  * nullptr or a pointer to a Configuration. This Configuration will be used
  * to retrieve the Solution from the inner Block of the BendersBFunction, at
  * every stage, after it is solved. This Configuration will be passed to
  * get_var_solution() of the inner Solver. The relevant part of the Solution
  * of the inner Block is the values of the active Variables of the
  * PolyhedralFunction. Thus, this Configuration can be used to specify that
  * only that portion of the Solution should be retrieved. Finally, the fourth
  * element, if present, must be either nullptr or a pointer to a
  * Configuration. This Configuration will be used to retrieve the dual
  * Solution from the inner Block of the BendersBFunction, at every stage,
  * after it is solved. This Configuration will be passed to
  * get_dual_solution() of the inner Solver.
  *
  * If the extra Configuration is not any of the specified above, an exception
  * is thrown.
  *
  * Here, we are assuming that the same Configuration can be applied to the
  * inner Block of the BendersBFunction at all stages. However, in principle,
  * the inner Block of the BendersBFunction at different stages could require
  * different Configuration. If this case ever happens, the implementation of
  * this method should be adapted to take it into consideration.
  *
  * If the given pointer to the ComputeConfig \p scfg is nullptr, then the
  * Configuration of this SDDPGreedySolver is reset to its default one.
  *
  * It is important to notice that every Configuration provided by \p scfg is
  * cloned (see Configuration::clone()) and, therefore, the caller is
  * responsible for destroying all these Configuration and the Configuration
  * pointed by \p scfg.
  *
  * @param scfg a pointer to a ComputeConfig.
  */

 void set_ComputeConfig( ComputeConfig *scfg = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the SDDPGreedySolver
 *  @{ */

 /// get the number of int parameters
 /** Get the number of int parameters.
  *
  * @return The number of int parameters.
  */

 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the number of string parameters
 /** Get the number of string parameters.
  *
  * @return The number of string parameters.
  */

 idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the number of vector-of-int parameters
 /** Get the number of vector-of-int  parameters.
  *
  * @return The number of vector-of-int parameters.
  */

 idx_type get_num_vint_par( void ) const override {
  return( idx_type( vintLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the number of vector-of-double parameters
 /** Get the number of vector-of-double  parameters.
  *
  * @return The number of vector-of-double parameters.
  */

 idx_type get_num_vdbl_par( void ) const override {
  return( idx_type( vdblLastAlgPar ) );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of an int parameter
 /** Get the default value of the int parameter with given index. Please see
  * the #int_par_type_SDDP_Greedy_S and #int_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 int get_dflt_int_par( const idx_type par ) const override {
  switch( par ) {
   case( intScenarioId ): return 0;
   case( intFirstStageScenarioId ): return -1;
   case( intUnregisterSolver ): return 0;
   case( intLogVerb ): return 0;
   case( intScenarioSeed ): return -1;
   case( intScenarioSampleFrequency ): return 1;
   case( intSimulationDataOutputPrecision ): return 20;
   case( intOutputScenario ): return -1;
   case( intEarlyConfig ): return 0;
   case( intLoadCutsOnce ): return 1;
  }
  return Solver::get_dflt_int_par( par );
 }

/*--------------------------------------------------------------------------*/

 /// get the default value of a string parameter
 /** Get the default value of the string parameter with given index. Please
  * see the #str_par_type_SDDP_Greedy_S and #str_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::string & get_dflt_str_par( const idx_type par ) const override {

  static const std::vector< std::string > default_values =
   { "" , "" , "" , "" , "" };

  if( par >= str_par_type_S::strLastAlgPar && par < strLastAlgPar )
   return default_values[ par - str_par_type_S::strLastAlgPar ];

  return Solver::get_dflt_str_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of a vector-of-int parameter
 /** Get the default value of the vector-of-int parameter with given index.
  * Please see the #vint_par_type_SDDP_Greedy_S and #vint_par_type_S
  * enumerations for a detailed explanation of the possible parameters. This
  * function returns the following values depending on the desired parameter:
  *
  * - #vintStagesSample: an empty vector
  *
  * For any other parameter, see Solver::get_dflt_vint_par().
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::vector< int > & get_dflt_vint_par( const idx_type par )
  const override {
  const static std::vector< int > empty;

  if( par == vintStagesSample ) {
   return empty;
  }

  return Solver::get_dflt_vint_par( par );
 }

/*--------------------------------------------------------------------------*/
 /// get the default value of a vector-of-double parameter
 /** Get the default value of the vector-of-double parameter with given index.
  * Please see the #vdbl_par_type_SDDP_Greedy_S and #vdbl_par_type_S
  * enumerations for a detailed explanation of the possible parameters. This
  * function returns the following values depending on the desired parameter:
  *
  * - #vdblInitialState: an empty vector
  *
  * For any other parameter, see Solver::get_dflt_vdbl_par().
  *
  * @param par The parameter whose default value is desired.
  *
  * @return The default value of the given parameter.
  */

 const std::vector< double > & get_dflt_vdbl_par( const idx_type par )
  const override {
  const static std::vector< double > empty;

  if( par == vdblInitialState ) {
   return empty;
  }

  return Solver::get_dflt_vdbl_par( par );
 }

/*--------------------------------------------------------------------------*/

 /// get a specific integer (int) numerical parameter
 /** Get a specific integer (int) numerical parameter. Please see the
  * #int_par_type_SDDP_Greedy_S and #int_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 int get_int_par( const idx_type par ) const override {
  switch( par ) {
   case( intScenarioId ): return f_scenario_id;
   case( intFirstStageScenarioId ): return f_first_stage_scenario_id;
   case( intUnregisterSolver ): return f_unregister_solver;
   case( intLogVerb ): return log_verbosity;
   case( intScenarioSeed ): return( f_seed == Inf< Index >() ) ? -1 : f_seed;
   case( intScenarioSampleFrequency ): return f_scenario_sample_frequency;
   case( intSimulationDataOutputPrecision ):
    return f_simulation_data_output_precision;
   case( intOutputScenario ): return f_output_scenario;
   case( intEarlyConfig ): return f_early_config;
   case( intLoadCutsOnce ): return f_load_cuts_once;
  }
  return( Solver::get_dflt_int_par( par ) );
 }

/*--------------------------------------------------------------------------*/

 /// get a specific string numerical parameter
 /** Get a specific string numerical parameter. Please see the
  * #str_par_type_SDDP_Greedy_S and #str_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::string & get_str_par( const idx_type par ) const override {
  switch( par ) {
   case( strInnerBC ): return f_inner_block_config_filename;
   case( strInnerBSC ): return f_inner_block_solver_config_filename;
   case( strLoadCuts ): return f_load_cuts_filename;
   case( strRandomCutsFile ): return f_random_cuts_filename;
   case( strSimulationData ): return f_simulation_data_filename;
  }
  return Solver::get_str_par( par );
 }

/*--------------------------------------------------------------------------*/

 /// get a specific vector-of-int parameter
 /** Get a specific vector-of-int parameter. Please see the
  * #vint_par_type_SDDP_Greedy_S and #vint_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::vector< int > & get_vint_par( const idx_type par )
  const override {
  switch( par ) {
   case( vintStagesSample ): return v_stages_to_sample;
  }
  return Solver::get_vint_par( par );
 }

/*--------------------------------------------------------------------------*/

 /// get a specific vector-of-double parameter
 /** Get a specific vector-of-double parameter. Please see the
  * #vdbl_par_type_SDDP_Greedy_S and #vdbl_par_type_S enumerations for a
  * detailed explanation of the possible parameters.
  *
  * @param par The parameter whose value is desired.
  *
  * @return The value of the given parameter.
  */

 const std::vector< double > & get_vdbl_par( const idx_type par )
  const override {
  switch( par ) {
   case( vdblInitialState ): return v_initial_state;
  }
  return Solver::get_vdbl_par( par );
 }

/*--------------------------------------------------------------------------*/

 /// returns the index of the int parameter with given string \p name
 /** This method takes a string, which is assumed to be the name of an int
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it. The method is given a void
  * implementation (throwing exception), rather than being pure virtual, so
  * that derived classes not having any int parameter do not have to bother
  * with implementing it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type int_par_str2idx( const std::string & name ) const override {
  if( name == "intScenarioId" ) return intScenarioId;
  if( name == "intFirstStageScenarioId" ) return intFirstStageScenarioId;
  if( name == "intUnregisterSolver" ) return intUnregisterSolver;
  if( name == "intScenarioSeed" ) return intScenarioSeed;
  if( name == "intScenarioSampleFrequency" ) return intScenarioSampleFrequency;
  if( name == "intSimulationDataOutputPrecision" )
   return intSimulationDataOutputPrecision;
  if( name == "intOutputScenario" ) return intOutputScenario;
  if( name == "intEarlyConfig" ) return intEarlyConfig;
  if( name == "intLoadCutsOnce" ) return intLoadCutsOnce;
  return Solver::int_par_str2idx( name );
 }

/*--------------------------------------------------------------------------*/

 /// returns the index of the string parameter with given string name
 /** This method takes a string, which is assumed to be the name of a string
  * parameter, and returns its index, i.e., the integer value that can be
  * used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type str_par_str2idx( const std::string & name ) const override {
  if( name == "strInnerBC" ) return strInnerBC;
  if( name == "strInnerBSC" ) return strInnerBSC;
  if( name == "strLoadCuts" ) return strLoadCuts;
  if( name == "strRandomCutsFile" ) return strRandomCutsFile;
  if( name == "strSimulationData" ) return strSimulationData;
  return Solver::str_par_str2idx( name );
 }

/*--------------------------------------------------------------------------*/

 /// returns the index of the vector-of-int parameter with given string name
 /** This method takes a string, which is assumed to be the name of a
  * vector-of-int parameter, and returns its index, i.e., the int value that
  * can be used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type vint_par_str2idx( const std::string & name ) const override {
  if( name == "vintStagesSample" ) return vintStagesSample;
  return Solver::vint_par_str2idx( name );
 }

/*--------------------------------------------------------------------------*/

 /// returns the index of the vector-of-double parameter with given string name
 /** This method takes a string, which is assumed to be the name of a
  * vector-of-double parameter, and returns its index, i.e., the double value
  * that can be used in [set/get]_par() to set/get it.
  *
  * @param name The name of the parameter.
  *
  * @return The index of the parameter with the given \p name.
  */

 idx_type vdbl_par_str2idx( const std::string & name ) const override {
  if( name == "vdblInitialState" ) return vdblInitialState;
  return Solver::vdbl_par_str2idx( name );
 }

/*--------------------------------------------------------------------------*/

/// returns the string name of the int parameter with given index
 /** This method takes an int parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & int_par_idx2str( const idx_type idx ) const override {

  static const std::vector< std::string > parameter_names =
   { "intScenarioId" , "intFirstStageScenarioId" , "intUnregisterSolver" ,
     "intScenarioSeed" , "intScenarioSampleFrequency" ,
     "intSimulationDataOutputPrecision" , "intOutputScenario" ,
     "intEarlyConfig" , "intLoadCutsOnce" };

  if( idx >= int_par_type_S::intLastAlgPar && idx < intLastAlgPar )
   return parameter_names[ idx - int_par_type_S::intLastAlgPar ];

  return Solver::int_par_idx2str( idx );
 }

/*--------------------------------------------------------------------------*/

 /// returns the string name of the string parameter with given index
 /** This method takes a string parameter index, i.e., the integer value that
  * can be used in [set/get]_par() [see above] to set/get it, and returns its
  * "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & str_par_idx2str( const idx_type idx ) const override {

  static const std::vector< std::string > parameter_names =
   { "strInnerBC" , "strInnerBSC" , "strLoadCuts" , "strRandomCutsFile" ,
     "strSimulationData" };

  if( idx >= str_par_type_S::strLastAlgPar && idx < strLastAlgPar )
   return parameter_names[ idx - str_par_type_S::strLastAlgPar ];

  return Solver::str_par_idx2str( idx );
 }

/*--------------------------------------------------------------------------*/

 /// returns the string name of the vector-of-double parameter with given index
 /** This method takes a vector-of-double parameter index, i.e., the double
  * value that can be used in [set/get]_par() [see above] to set/get it, and
  * returns its "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & vint_par_idx2str( const idx_type idx ) const override {
  static const std::vector< std::string > parameter_names =
   { "vintStagesSample" };
  if( idx >= vint_par_type_S::vintLastAlgPar && idx < vintLastAlgPar )
   return parameter_names[ idx - vint_par_type_S::vintLastAlgPar ];
  return Solver::vint_par_idx2str( idx );
 }

/*--------------------------------------------------------------------------*/

 /// returns the string name of the vector-of-double parameter with given index
 /** This method takes a vector-of-double parameter index, i.e., the double
  * value that can be used in [set/get]_par() [see above] to set/get it, and
  * returns its "string name".
  *
  * @param idx The index of the parameter.
  *
  * @return The name of the parameter with the given index \p idx.
  */

 const std::string & vdbl_par_idx2str( const idx_type idx ) const override {
  static const std::vector< std::string > parameter_names =
   { "vdblInitialState" };
  if( idx >= vdbl_par_type_S::vdblLastAlgPar && idx < vdblLastAlgPar )
   return parameter_names[ idx - vdbl_par_type_S::vdblLastAlgPar ];
  return Solver::vdbl_par_idx2str( idx );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the model encoded by the current Block
 *  @{ */

 /// (try to) solve the model encoded in the SDDPBlock for a single scenario
 /** This method tries to solve the deterministic (single-scenario) multistage
  * problem defined in (2) above. The problem is determined by the scenario
  * whose id is returned by the get_scenario_id() method and whose initial
  * state is defined in the SDDPBlock. This method does not really try to
  * solve the problem (2), but employs a procedure that may find a feasible
  * solution for that problem and can be interpreted as a simulation.
  *
  * Beginning at the first stage, this method tries to solve the subproblem
  * associated with each stage, sequentially, until the last one. The
  * subproblem at stage t is encoded by the inner Block of the
  * BendersBFunction associated with the stage t. The subproblem at stage t >
  * 0 may depend on the variables of the subproblem at stage t-1. After
  * solving the subproblem at stage t-1, the solution found to this subproblem
  * is used to update the next subproblem according to this dependency (which
  * is characterized by the BendersBFunction).
  *
  * At any given stage, the subproblem may be successfully solved or not. If a
  * solution to a subproblem is not found (for instance, if the subproblem
  * turns out to be infeasible, unbounded, or an error occurs while solving
  * it), then this method stops with the corresponding status as described in
  * #sddp_greedy_sol_type. In this case, no solution for the deterministic
  * (single-scenario) multistage problem can be provided.
  *
  * Notice that a feasible solution may not be found for the deterministic
  * (single-scenario) multistage problem (2) even if one exists.
  *
  * @return Please refer to #sddp_greedy_sol_type for a description of each
  *         value that this method may return.
  */

 int compute( bool changedvars = true ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR CHANGING THE DATA OF THE SDDPGreedySolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the SDDPGreedySolver
 *  @{ */

 /// sets the scenario that should be considered
 /** This method defines which scenario should be considered when trying to
  * solve the deterministic single-scenario problem. The \p scenario_id
  * parameter must be the id of a scenario handled by the SDDPBlock.
  *
  * @param scenario_id The id of the scenario. */

 void set_scenario_id( Index scenario_id ) {
  if( f_scenario_id == scenario_id )
   return;
  f_scenario_id = scenario_id;
  status_compute = Solver::kUnEval;
 }

/*--------------------------------------------------------------------------*/

 /// sets the scenario that should be considered at the first stage
 /** This method defines which scenario should be considered at the first
  * stage when trying to solve the deterministic single-scenario problem. The
  * \p scenario_id parameter must be one of the following:
  *
  * - the id of a scenario handled by the SDDPBlock, which will then be the
  *   scenario considered at the first stage;
  *
  * - -1, which means that the id of the scenario considered at the first
  *   stage is that given by get_scenario_id();
  *
  * - any other negative number, which means that no scenario must be set
  *   while solving the sub-problem at the first stage (i.e., the data for the
  *   first stage subproblem has already been set).
  *
  * @param scenario_id The id of a scenario or a negative number. */

 void set_first_stage_scenario_id( int scenario_id ) {
  if( f_first_stage_scenario_id == scenario_id )
   return;
  f_first_stage_scenario_id = scenario_id;
  status_compute = Solver::kUnEval;
 }

/*--------------------------------------------------------------------------*/

 /// sets the scenario to be considered
 /** This method updates the sub-Block of the SDDPBlock associated with the
  *  given \p stage with the data provided by the scenario whose id is
  *  \p scenario_id.
  *
  * @param scenario_id The id of a scenario handled by the SDDPBlock.
  *
  * @param stage An integer between 0 and get_time_horizon() - 1.
  */
 void set_scenario( Index scenario_id , Index stage );

/*--------------------------------------------------------------------------*/

 /// sets the callback function
 /** It sets the callback function that is called right before the sub-problem
  * at each stage is solved. The parameter of the callback function is the
  * stage associated with the sub-problem that will be solved.
  */
 void set_callback( std::function< void( Index ) > function ) {
  callback = function;
 }

/*--------------------------------------------------------------------------*/

 /// sets the random number engine used to select the scenarios
 /** This function sets the random number engine that is used to select the
  * scenario at each stage when random scenarios must be considered (see the
  * #intScenarioSeed parameter).
  *
  * @param A random number engine.
  */
 void set_random_number_engine( std::mt19937 random_number_engine ) {
  this->random_number_engine = random_number_engine;
 }

/**@} ----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any)
 * @{ */

 bool has_var_solution( void ) override {
  return f_has_var_solution;
 }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool has_dual_solution( void ) override;

/*--------------------------------------------------------------------------*/

 bool is_dual_feasible( void ) override;

/*--------------------------------------------------------------------------*/

 void get_dual_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_dual_solution( void ) override;

/*--------------------------------------------------------------------------*/

 bool has_dual_direction( void ) override;

/*--------------------------------------------------------------------------*/

 void get_dual_direction( Configuration * dirc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 bool new_dual_direction( void ) override;

/*--------------------------------------------------------------------------*/

 OFValue get_var_value( void ) override {
  return solution_value;
 }

/*--------------------------------------------------------------------------*/

 OFValue get_lb( void ) override {
  if( ( get_objective_sense() == Objective::eMax ) && has_var_solution() )
   return solution_value;
  return( - Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

 OFValue get_ub( void ) override {
  if( ( get_objective_sense() == Objective::eMin ) && has_var_solution() )
   return solution_value;
  return( Inf< OFValue >() );
 }

/*--------------------------------------------------------------------------*/

 /** If a call to compute() returns kError, kStopTime, or kStopIter, this
  * method returns the stage at which the associated event has
  * occurred. Otherwise, this method returns Inf< Index >().
  *
  * @return the stage at which a fault has occurred.
  */
 Index get_fault_stage() const {
  return fault_stage;
 }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE SDDPGreedySolver -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the state of the SDDPGreedySolver
 *  @{ */

 /// returns the time horizon of the SDDPBlock attached to this SDDPGreedySolver
 /** Returns the time horizon of the problem represented by the SDDPBlock
  * attached to this SDDPGreedySolver.
  *
  * @return The time horizon of the problem represented by the SDDPBlock
  *         attached to this SDDPGreedySolver.
  */
 Index get_time_horizon( void ) const {
  if( f_Block )
   return static_cast< SDDPBlock * >( f_Block )->get_time_horizon();
  return 0;
 }

/*--------------------------------------------------------------------------*/

 /// returns the id of the scenario being currently considered
 /** This method returns the id of the scenario being currently considered
  * when trying to solve the deterministic (single-scenario) multistage
  * problem.
  *
  * @return The id of the scenario to be considered.
  */
 Index get_scenario_id( void ) const {
  return f_scenario_id;
 }

/*--------------------------------------------------------------------------*/

 /// returns the id of the scenario considered at the first stage
 /** This method returns the id of the scenario considered at the first stage
  * when trying to solve the deterministic (single-scenario) multistage
  * problem or a negative number as follows:
  *
  * - If the value returned is nonnegative, then it is the id of a scenario
  *   handled by the SDDPBlock, which is the scenario considered at the first
  *   stage.
  *
  * - If the value returned is -1, then scenario considered at the first stage
  *   is that given by get_scenario_id().
  *
  * - If the value returned is any other negative number then this means that
  *   no scenario must be set while solving the sub-problem at the first
  *   stage.
  *
  * @return The id of the scenario to be considered or a negative number.
  */
 Index get_first_stage_scenario_id( void ) const {
  return f_first_stage_scenario_id;
 }

/*--------------------------------------------------------------------------*/

 /// returns the status of the most recent call to compute()
 /** Returns the status of the most recent call to compute().
  *
  * @return The status of the most recent call to compute().
  */
 Index get_status( void ) const {
  return status_compute;
 }

/*--------------------------------------------------------------------------*/

 /// returns a copy of the random number engine used to select the scenarios
 /** This function returns a copy of the random number engine that is used to
  * select the scenario at each stage when random scenarios must be considered
  * (see the #intScenarioSeed parameter).
  *
  * @return A copy of the random number engine used to select the scenarios.
  */
 std::mt19937 get_random_number_engine() const {
  return random_number_engine;
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 /// The id of the scenario that should be considered
 /**< This is the id of the scenario that must be considered at all stages
  * (except possibly the first stage). */
 Index f_scenario_id = 0;

 /// The id of the scenario that must be considered at the first stage
 int f_first_stage_scenario_id = -1;

 /// It determines how (and if) the scenario should be output
 int f_output_scenario = -1;

 /// The stage at which some special event has happened
 Index fault_stage = Inf< Index >();

 /// It indicates whether the inner Blocks should be configured in advance
 int f_early_config = 0;

 /// It indicates whether cuts should be loaded only once (see #intLoadCutsOnce)
 int f_load_cuts_once = 1;

 /// Indicates whether the Solver of the inner Block must be unregister
 bool f_unregister_solver = false;

 /// Function to be called right before each sub-problem is solved
 std::function< void( Index ) > callback;

 /// Name of the default BlockConfig file for the inner Blocks
 std::string f_inner_block_config_filename{};

 /// Default BlockConfig for the inner Blocks
 BlockConfig * f_inner_block_config = nullptr;

 /// Name of the default BlockSolverConfig file for the inner Blocks
 std::string f_inner_block_solver_config_filename{};

 /// Default BlockConfig for the inner Blocks
 BlockSolverConfig * f_inner_block_solver_config = nullptr;

 /// Configuration to be passed to the get_var_solution() method
 Configuration * f_get_var_solution_config = nullptr;

 /// Configuration to be passed to the get_dual_solution() method
 Configuration * f_get_dual_solution_config = nullptr;

 /// Names of the BlockConfig file for the inner Blocks
 std::vector< std::string > v_BC_filename;

 /// BlockConfig for the inner Blocks
 std::vector< BlockConfig * > v_BC;

 /// Names of the BlockSolverConfig file for the inner Blocks
 std::vector< std::string > v_BSC_filename;

 /// BlockSolverConfig for the inner Blocks
 std::vector< BlockSolverConfig * > v_BSC;

 /// Indicates whether the inner Block of each BendersBFunction was configured
 std::vector< bool > v_inner_block_configured;

 /// Indicates whether the Solver of the inner Block of each BendersBFunction
 /// has been configured
 std::vector< bool > v_inner_solver_configured;

 /// It indicates whether cuts for each stage have been loaded
 std::vector< bool > v_cuts_loaded;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE CLASSES -------------------------------*/
/*--------------------------------------------------------------------------*/

 class SimulationData {

  using matrix = std::vector< std::pair< Index , std::vector< double > > >;

 public:

  void store_initial_state( const std::vector< double > & initial_state ,
                            Index stage ) {
   initial_states.emplace_back( stage , initial_state );
  }

  void store_subgradient_final_state( std::vector< double > && subgradient ,
                                      Index stage ) {
   subgradients_final_state.emplace_back( stage , std::move( subgradient ) );
  }

  void store_subgradient_initial_state( std::vector< double > && subgradient ,
                                        Index stage ) {
   subgradients_initial_state.emplace_back( stage , std::move( subgradient ) );
  }

  void store_objective_value( double objective_value , Index stage ) {
   objective_values.emplace_back( stage , objective_value );
  }

  void clear() {
   subgradients_initial_state.clear();
   subgradients_final_state.clear();
   initial_states.clear();
   objective_values.clear();
  }

  // Subgradients of the objective of the subproblems with respect to the
  // initial state
  matrix subgradients_initial_state;

  // Subgradients of the objective of the subproblems with respect to the
  // final state
  matrix subgradients_final_state;

  // Initial states
  matrix initial_states;

  // Values of the objectives of the subproblems disregarding the future value
  // function
  std::vector< std::pair< Index , double > > objective_values;
 };

/*--------------------------------------------------------------------------*/

 class Logger {

 public:

  Logger( SDDPGreedySolver * solver ,
          std::ostream * log_stream , int log_verbosity ) {
   this->solver = solver;
   this->log_verbosity = log_verbosity;
   this->f_log = log_stream;

   if( solver ) {
    this->stage_width =
     std::max( std::string::size_type( 6 ) ,
               std::to_string( solver->get_time_horizon() ).size() );
    const auto sddp_block = static_cast< SDDPBlock * >( solver->get_Block() );
    this->scenario_width =
     std::max( std::string::size_type( 9 ) ,
               std::to_string( sddp_block->get_scenario_set().size() ).size() );
   }
  }

  void log( double objective_value , double future_value ) const;
  void log( double objective_value ) const;
  void log() const;
  void log( Index stage , Index scenario ) const;
  void log_header() const;
  void show_status() const;

 private:

  const int precision = 7;
  const int width = precision + 6;
  unsigned long stage_width;
  unsigned long scenario_width;
  int log_verbosity = 0;
  std::ostream * f_log = nullptr;
  SDDPGreedySolver * solver = nullptr;
 };

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// solves the subproblem associated with the given stage
 /** This method solves the subproblem associated with the given \p stage. If
  * \p write_solution is true then the solution found for the subproblem (if
  * any) is written into its Block. To solve the subproblem, the method
  * compute() of the Solver attached to its Block is invoked and the status
  * returned by that method is returned here.
  *
  * @param stage The stage associated with the subproblem to be solved.
  *
  * @param write_solution Indicates whether the solution found for the
  *        subproblem (if any) must be written into its Block.
  *
  * @return The status returned by compute() when solving the subproblem.
  */
 int solve( Index stage , bool write_solution = false );

/*--------------------------------------------------------------------------*/

 /// returns the Solver attached to the subproblem at the given stage
 /** This method returns a pointer to the Solver attached to the subproblem
  * associated with the given \p stage.
  *
  * @param stage The stage associated with the subproblem whose Solver is
  *        desired.
  *
  * @return A pointer to the Solver attached to the subproblem associated with
  *         the given \p stage.
  */
 Solver * get_sub_solver( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// returns a pointer to the BendersBFunction associated with the given stage
 /** This method returns a pointer to the BendersBFunction associated with the
  * given \p stage.
  *
  * @param[in] stage An Index in the interval [0, T-1], where T is the time
  *            horizon.
  *
  * @return A pointer to the BendersBFunction associated with the given stage.
  */
 BendersBFunction * get_benders_function( Index stage ) const {
  if( stage >= get_time_horizon() )
   throw( std::invalid_argument( "SDDPGreedySolver::get_benders_function: "
                                 "invalid stage index: " +
                                 std::to_string( stage ) ) );

  auto benders_block = static_cast< BendersBlock * >
   ( static_cast< SDDPBlock * >( f_Block )->
     get_sub_Block( stage )->get_inner_block() );

  auto objective = static_cast< FRealObjective * >
   ( benders_block->get_objective() );

  return static_cast< BendersBFunction * >( objective->get_function() );
 }

/*--------------------------------------------------------------------------*/

 /// returns the solution associated with the problem at the given stage
 /** This function returns the solution of the problem associated with the
  * given \p stage, which is part of the state variables of the next stage.
  *
  * @param stage The stage whose solution is required.
  *
  * @return The vector containing the solution of the problem at the given
  *         stage.
  */
 std::vector< double > get_solution( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// sets the state variables of the subproblem at the given stage
 /** This function sets the state variables associated with the subproblem at
  * the given \p stage.
  *
  * @param The vector containing the state of the problem at the given stage.
  *
  * @param stage The stage whose state must be set.
  */
 void set_state( const std::vector< double > & state , Index stage ) const;

/*--------------------------------------------------------------------------*/

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/

 double get_sub_solution_value( Index stage ) const {
  assert( stage < get_time_horizon() );
  double solution_value = 0;
  auto sub_solver = get_sub_solver( stage );
  if( sub_solver->is_var_feasible() )
   solution_value = sub_solver->get_var_value();
  else if( get_objective_sense( stage ) == Objective::eMin )
   solution_value = sub_solver->get_ub();
  else
   solution_value = sub_solver->get_lb();
  return solution_value - get_future_value( stage );
 }

/*--------------------------------------------------------------------------*/

 double get_future_value( Index stage ) const {
  assert( stage < get_time_horizon() );
  return static_cast< SDDPBlock * >( f_Block )->get_future_cost( stage , 0 );
 }

/*--------------------------------------------------------------------------*/

 int get_objective_sense( Index stage = 0 ) const {
  assert( stage < get_time_horizon() );
  auto benders_function = get_benders_function( stage );
  auto inner_block = benders_function->get_inner_block();
  assert( inner_block );
  return inner_block->get_objective_sense();
 }

/*--------------------------------------------------------------------------*/

 void configure_inner_block( Index stage );

/*--------------------------------------------------------------------------*/

 void unregister_solver_inner_block( Index stage );

/*--------------------------------------------------------------------------*/

 /// output the simulation data (if any)
 /** This function outputs the data obtained during the simulation (i.e., the
  * most recent call to compuet()), which includes subgradients with respect
  * to the initial and final states, objective values, initial states, and
  * scenarios, into the file with the given \p filename. The format of the
  * file will follow that specified in the description of #strSimulationData.
  *
  * @param filename The name of the file in which the data obtained during the
  *        simulation should be stored. */

 void output_simulation_data( const std::string & filename ) const;

/*--------------------------------------------------------------------------*/

 /// load cuts for the given \p stage
 /** This function loads cuts for the given \p stage from a file that can be
  * specified by the #strLoadCuts parameter. If no file has been specified,
  * this function does nothing. If a file has been specified, cuts for the
  * given \p stage will be retrieved from that file and loaded in the
  * corresponding PolyhedralFunction(s). If the SDDPBlock has multiple
  * sub-Blocks per stage then the cuts are loaded in every sub-Block. However,
  * this function assumes that there is only one PolyhedralFunction per
  * sub-Block. See the #strLoadCuts parameter for a description of the format
  * that the file must have.
  *
  * It is important to notice that the cuts are added to the
  * PolyhedralFunction and any other cuts that were possibly already there in
  * the PolyhedralFunction are kept there. Moreover, cuts for stage t are
  * loaded within compute(), right before the subproblem associated with stage
  * t is solved. If cuts are not to be loaded on subsequent calls to
  * compute(), the value of parameter #strLoadCuts must be updated to the
  * empty string.
  *
  * @param stage The stage (between 0 and get_time_horizon() - 1) for which
  *        cuts should be loaded. */

 void load_cuts( Index stage );

/*--------------------------------------------------------------------------*/

 /// returns the id of the scenario that must be considered at the given stage
 /** This function returns the id of the scenario that must be considered at
  * the given \p stage.
  *
  * @param stage A stage (between 0 and get_time_horizon() - 1).
  *
  * @return The id of the scenario that must be considered at the given
  *         \p stage. */

 Index get_scenario_id( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// checks whether a scenario should be sampled for the given \p stage
 /** This function returns true if and only if a scenario should be sampled
  * for the given \p stage.
  *
  * @param stage A stage between 0 and get_time_horizon() - 1.
  *
  * @return true if and only if a scenario should be sampled for the given \p
  *         stage.
  */
 bool should_sample( Index stage ) const;

/*--------------------------------------------------------------------------*/

 /// sample a new scenario for the given \p stage
 /** This function selects, at random, a new scenario for the given \p stage.
  *
  * @param stage A stage (between 0 and get_time_horizon() - 1). */

 void sample_scenario( Index stage );

/*--------------------------------------------------------------------------*/

 void store_simulation_data( Index stage , Index scenario_index );

/*--------------------------------------------------------------------------*/

 void store_subgradient_final_state( Index stage );

/*--------------------------------------------------------------------------*/

 void store_subgradient_initial_state
  ( Index stage , Index scenario_index ,
    const std::vector< double > & initial_state );

/*--------------------------------------------------------------------------*/

 void reset_compute_time() {
  f_compute_start_time = std::chrono::system_clock::now();
 }

/*--------------------------------------------------------------------------*/

 double get_compute_time() const {
  const auto now = std::chrono::system_clock::now();
  std::chrono::duration< double > time = now - f_compute_start_time;
  return time.count();
 }

/*--------------------------------------------------------------------------*/

 void reset_subproblem_time() {
  f_subproblem_start_time = std::chrono::system_clock::now();
 }

/*--------------------------------------------------------------------------*/

 double get_subproblem_time() const {
  const auto now = std::chrono::system_clock::now();
  std::chrono::duration< double > time = now - f_subproblem_start_time;
  return time.count();
 }

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// The status returned by compute()
 int status_compute = Solver::kUnEval;

 /// It indicates whether a feasible solution has been found
 bool f_has_var_solution = false;

 /// It indicates whether all of the sub-Solvers have dual solutions
 bool f_has_dual_solution = false;

 /// The value of the solution (if any)
 double solution_value = 0.0;

 /// It indicates the level of verbosity of the log
 int log_verbosity = 0;

 /// Precision of the simulation data that is output
 int f_simulation_data_output_precision = 20;

 /// The name of the file out of which cuts are loaded
 std::string f_load_cuts_filename;

 /// The name of the file containing the random cuts
 std::string f_random_cuts_filename;

 /// The name of the file in which the simulation data will be saved
 std::string f_simulation_data_filename;

 /// Initial state for the first stage problem
 std::vector< double > v_initial_state;

 /// Data obtained during the simulation
 SimulationData f_simulation_data;

 /// Frequency at which scenarios should be sampled
 int f_scenario_sample_frequency = 1;

 /// IDs of the scenarios to be considered at each stage
 std::vector< Index > v_random_scenario_id;

 /// Stages at which scenarios must be sampled
 std::vector< int > v_stages_to_sample;

 /// Distribution for selecting the scenarios at each stage
 std::uniform_int_distribution< Index > scenario_distribution;

 /// Random number engine to select the scenarios
 std::mt19937 random_number_engine;

 /// Seed for the random number engine that selects the scenarios
 Index f_seed = Inf< Index >();

 /// Point in time at which the most recent call to compute() has started
 std::chrono::time_point< std::chrono::system_clock > f_compute_start_time;

 /// Point in time at which the solution of the most recent subproblem started
 std::chrono::time_point< std::chrono::system_clock > f_subproblem_start_time;

};   // end( class SDDPGreedySolver )

/** @} end( group( SDDPGreedySolver_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* SDDPGreedySolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File SDDPGreedySolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
