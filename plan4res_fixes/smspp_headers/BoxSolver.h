/*--------------------------------------------------------------------------*/
/*-------------------------- File BoxSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BoxSolver class, which implements a CDASolver for
 * problems (or relaxations thereof) with an extremely simple structure:
 * only bound (box) Constraint on the ColVariable and a separable Objective
 * (a FRealObjective with either a LinearFunction or a DQuadFunction inside).
 * This only rarely is a significant problem in itself (although it may
 * indeed appear in the context of relaxation methods where everything else
 * has been relaxed), but it can be used to quickly obtain (hopefully,
 * finite) bounds on the optimal value of more complex problems that may be
 * useful during algorithmic approaches. For this reason, BoxSolver has an
 * uncommon "lax" attitude w.r.t. all the Constraint in the Block that are
 * not box ones: rather than protesting for their existance and refusing to
 * load, it plainly ignores them. This means that the computed optimal value
 * possibly is a(n hopefully finite) valid bound (lower or upper, according
 * to the sense of the original Objective) on the true optimal value. Yet,
 * this bound is obtained quickly.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BoxSolver
 #define __BoxSolver  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class Function;          // forward declaration of Function

 class OneVarConstraint;  // forward declaration of OneVarConstraint

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BoxSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Solver for "extremely simple" Block, or relaxations thereof
/** The BoxSolver class implements the CDASolver interface for problems (or
 * relaxations of problems) with an extremely simple structure:
 *
 * - only ColVariable
 *
 * - only bound (box) Constraint on the ColVariable, i.e., either those that
 *   are "inherent" in the ColVariable or :OneVarConstraint;
 *
 * - a separable Objective, i.e., a FRealObjective with either a
 *   LinearFunction or a DQuadFunction inside.
 *
 * This only rarely is a significant problem in itself (although it may
 * indeed appear in the context of relaxation methods where everything else
 * has been relaxed), but it can be used to quickly obtain (hopefully,
 * finite) bounds on the optimal value of more complex problems that may be
 * useful during algorithmic approaches. For this reason, BoxSolver has an
 * uncommon "lax" attitude w.r.t. all the Constraint in the Block that are
 * not box ones: rather than protesting for their existance and refusing to
 * load, it plainly ignores them. This means that the computed optimal value
 * possibly is a(n hopefully finite) valid bound (lower or upper, according
 * to the sense of the original Objective) on the true optimal value. Yet,
 * this bound is obtained quickly.
 *
 * In fact, for the problems solved by BoxSolver have the very uncommon
 * property that it is basically as costly to compute either the min or the
 * max of the Objective, and to compute *both* the min and the max. So, this
 * is what BoxSolver does: each time compute() is called, it computes and
 * makes it available both the minimum and the maximum of the objective,
 * (almost) regardless to what the original sense was.
 *
 * Note that, albeit very simple due to separability, the problem is not
 * necessarily convex: the ColVariable may have integrality constraints, and
 * a quadratic function is not equivalent to minimize or maximise as it is a
 * linear one. Thus, while BoxSolver derives from CDASolver,
 *
 *    IT WILL ONLY REPORT DUAL SOLUTIONS FOR THOSE ColVariable WHOSE
 *    OPTIMIZATION LEADS TO EXACT DUALITY (no integrality constraints,
 *    minimization of convex / maximization of concave.
 *
 * Furthermore
 *
 *    DUAL VALUES CAN ONLY BE RETURNED IF THE "ACTIVE" BOX CONSTRAINT IS
 *    A OneVarConstraint, WHICH MAY NOT ALWAYS BE SO
 *
 * In fact, ColVariable has "built-in" bound constraints which may be the
 * ones dictating the optimal solution (the box constraints dictated by the
 * OneVarConstraint being looser, ar not being there at all). In this case,
 * the optimal dual value "has nowhere to go" and is lost for good. If this
 * is not acceptable, the user must ensure that a OneVarConstraint always
 * exists that is at least as tight as the "built-in" bound constraints. */

class BoxSolver : public CDASolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */

 /// public enum "extending" int_par_type_CDAS to BoxSolver

 enum int_par_type_BoxS {
  intPDSol = intLastParCDAS ,   ///< what is computed besides the bound

  intLastParBoxS    ///< first allowed parameter value for derived classes
                    /**< convenience value for easily allow derived classes
                     * to further extend the set of types of return codes */
  };

/*--------------------------------------------------------------------------*/
/*----------------- CONSTRUCTING AND DESTRUCTING BoxSolver -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing BoxSolver
 *  @{ */

 /// empty constructor

 BoxSolver( void ) : CDASolver() , f_sol( 0 ) , f_sol_comp( 0 ) ,
  f_state( kUnEval ) , f_max_val( - Inf< OFValue >() ) ,
  f_min_val( Inf< OFValue >() ) , f_feas( false ) , f_sense( -1 ) ,
  f_altobj( nullptr ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: it really does nothing since v_mod is empty

 virtual ~BoxSolver() {}

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set an alternative Function as the (unique) Objective of the problem
 /** This method instructs BoxSolver to temporarily ignore all the Objective
  * in the Block and rather do its computation using the provided Function as
  * if it were the only Objective. Currently, the only accepted Function are
  *
  * - LinearFunction
  *
  * - DQuadFunction
  *
  * Note that since all the other Objective are disregarded, each Variable in
  * the Block that is not active un \p newf does not contribute to the
  * computation of the optimal values and solutions.
  *
  * By passing nullptr the BoxSolver is instructed to return considering the
  * Objective in the Block (and all its sub-Block, recursively). */

 void set_Objective_Function( Function * newf ) {
  if( newf != f_altobj ) {
   f_altobj = newf;
   f_state = kUnEval;
   }
  }

/*--------------------------------------------------------------------------*/
 /// set the (unique) integer parameter of BoxSolver, see set_sol()

 /* FIX (overload hiding): declaring any set_par() here hides ALL base-class
  * set_par() overloads. With no double overload and no using-declaration, a
  * set_par( par , <double> ) call narrows the double to int and dispatches to
  * the int overload - whose enum is a SEPARATE index space - silently
  * corrupting an unrelated parameter. Same defect as MILPSolver, where
  * dblAAccSol (double idx 7) overwrote intRelaxIntVars (int idx 7). */
 using CDASolver::set_par;
 void set_par( idx_type par , int value ) override {
  if( par == intPDSol ) {
   set_sol( char( value ) );
   return;
   }
  CDASolver::set_par( par , value );
  }

/*--------------------------------------------------------------------------*/
 /// set what is computed (primal/dual solution)
 /** Decides what is computed besides the optimal value: the primal and/or
  * dual solution.
  *
  * Since computing solutions has the same cost as producing the bound(s), it
  * is better done inside compute(); in this case the optimal primal and / or
  * dual solution is immediately written in the ColVariable / dual values of
  * the :OneVarConstraint.
  *
  * @param sol is a char, coded bit-wise, that decides what is computed:
  *
  *        - bit 0: the primal solution is produced
  *
  *        - bit 1: the dual solution of the OneVarConstraint is produced
  *          (this is the only "real" part, since only the bounds are
  *          really taken into account during optimization)
  *
  *        - bit 2: the dual solution of all the other FRowConstraint is
  *          "produced", in the sense that it is set to 0 (to reflect the
  *          fact that the other constraint have been relaxed)
  *
  * Note that there is the usual issue with dual solutions: if one of the
  * "inherent" bounds of a ColVariable has a nonzero dual value but there is
  * no :OneVarConstraint with the same bound, then there is no place where to
  * store the dual value and it is "lost". If this is a problem, the Block
  * must always have the bounds specified via the :OneVarConstraint.
  *
  * However, doing this is a "minor" breach of the CDASolver interface, in
  * that primal / dual solutions should ony be written when
  * get_var_solution() and get_dual_solution() are called. This is why the
  * default is 0. If the primal / dual solution has not been written during
  * compute() it can still be required "normally" by calling
  * get_var_solution() / get_dual_solution() as per the normal interface.
  * This has roughly the same cost as solving the problem in the first place,
  * but there you go (anyway the cost cannot reasonably be painted as
  * "large"). */

 void set_sol( char sol = 0 ) { f_sol = sol & 7; }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the model encoded by the current Block
 *  @{ */

 /// solve the model encoded in the Block
 /** BoxSolver not only solves the model encoded in the Block (or a
  * relaxation thereof if there are other Constraint save the "box" ones),
  * but at the same time computes the maximum and the minimum of the
  * Objective over the feasible region. Of course the two can be infinite.
  * Also, compute() immediately writes the primal/dual optimal solution in
  * the ColVariable / dual values of the :OneVarConstraint, if so
  * instructed [see set_sol()]. */

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] OFValue get_lb( void ) override {
  if( f_state == kUnEval )
   throw( std::logic_error( "BoxSolver: compute() not called" ) );

  return( f_sense == 1 ? f_max_val : f_min_val );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] OFValue get_ub( void ) override {
  if( f_state == kUnEval )
   throw( std::logic_error( "BoxSolver: compute() not called" ) );

  return( f_sense == 1 ? f_max_val : f_min_val );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] bool has_var_solution( void ) override {
  return( ( f_state == kOK ) || ( f_state == kUnbounded ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] bool has_dual_solution( void ) override {
  return( ( f_state == kOK ) || ( f_state == kInfeasible ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] OFValue get_var_value( void ) override {
  if( f_state == kUnEval )
   throw( std::logic_error( "BoxSolver: compute() not called" ) );

  return( f_sense == 1 ? f_max_val : f_min_val );
  }

/*--------------------------------------------------------------------------*/

 void get_var_solution( Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_dual_solution( Configuration *solc = nullptr ) override;

 /*--------------------------------------------------------------------------*/

 [[nodiscard]] bool has_var_direction( void ) override {
  return( f_state == kUnbounded );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 /// one day perhaps we will implement it

 [[nodiscard]] bool has_dual_direction( void ) override {
  return( f_state == kInfeasible );
  }

----------------------------------------------------------------------------*/

 void get_var_direction( Configuration * dirc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ---
 /// one day perhaps we will implement it

 void get_dual_direction( Configuration * dirc = nullptr ) override;

----------------------------------------------------------------------------*/
 /// return the "opposite" bound w.r.t. get_var_value()
 /** While get_var_value() returns the value corresponding to the optimization
  * with the sense specified by the Objective, get_opposite_value() returns
  * the value corresponding to the optimization with the opposite sense. */

 OFValue get_opposite_value( void ) const {
  if( f_state == kUnEval )
   throw( std::logic_error( "BoxSolver: compute() not called" ) );

  return( f_sense == 0 ? f_max_val : f_min_val );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the Solver
 *  @{ */

 /// returns the current value controlling what is computed

 [[nodiscard]] char get_sol( void ) const { return( f_sol ); }

/*--------------------------------------------------------------------------*/


 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( CDASolver::get_num_int_par() + 1 );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  return( par == intPDSol ? 0 : CDASolver::get_dflt_int_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  return( par == intPDSol ? f_sol : CDASolver::get_int_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  return( name == "intPDSol" ? intPDSol
	                     : CDASolver::int_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static std::string _parnm( "intPDSol" );
  return( idx == intPDSol ? _parnm : CDASolver::int_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the model
 *  @{ */

 void add_Modification( sp_Mod & mod ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void check_sense( Block * blck );

/*--------------------------------------------------------------------------*/

 void reset( void ) { f_state = kUnEval; }

 void resetf( void ) { reset(); f_feas = false; }

/*--------------------------------------------------------------------------*/
 // returns true if var belongs to f_Block

 bool is_mine( ColVariable & var ) {
  for( auto bk = var.get_Block() ; bk ; bk = bk->get_f_Block() )
   if( bk == f_Block )
    return( true );

  return( false );
  }

/*--------------------------------------------------------------------------*/

 using VarValue = ColVariable::VarValue;

 /// get the bounds l and u of var
 void get_var_data( ColVariable & var , VarValue & l , VarValue & u );

 /// get the bounds l and u and the coefficients a and b of var
 void get_var_data( ColVariable & var , VarValue & l , VarValue & u ,
		    OFValue & a , OFValue & b );

 /// get the bounds l and u and their pointers of var
 void get_var_data( ColVariable & var , VarValue & l , VarValue & u ,
		    OneVarConstraint * & cl , OneVarConstraint * & cu );

 /// get the bounds l and u, their pointers and the coefficients a and b of var
 void get_var_data( ColVariable & var , VarValue & l , VarValue & u ,
		    OFValue & a , OFValue & b ,
		    OneVarConstraint * & cl , OneVarConstraint * & cu );

 /// solve a "linear" variable
 void sol_variable( ColVariable & var , VarValue l , VarValue u ,
		    OFValue b );

 /// solve a "quadratic" variable
 void sol_variable( ColVariable & var , VarValue l , VarValue u ,
		    OFValue a , OFValue b );

 /// solve a "linear" variable with bounds pointers
 void sol_variable( ColVariable & var , VarValue l , VarValue u ,
		    OFValue b ,
		    OneVarConstraint * cl , OneVarConstraint * cu );

 /// solve a "quadratic" variable with bounds pointers
 void sol_variable( ColVariable & var , VarValue l , VarValue u ,
		    OFValue a , OFValue b ,
		    OneVarConstraint * cl , OneVarConstraint * cu );

 /// just set any finite value to var, if any
 void process_variable_bnds( ColVariable & var );

 /// solve var getting all of its data
 void process_variable( ColVariable & var );

 /// produce a primal scolution for var getting all of its data
 void process_variable_sol( ColVariable & var );

 /// produce a primal solution to a  "linear" variable
 void process_var_sol( ColVariable & var , VarValue l , VarValue u ,
		       OFValue b );

 /// produce a primal solution to a  "quadratic" variable
 void process_var_sol( ColVariable & var , VarValue l , VarValue u ,
		       OFValue a , OFValue b );

 /// produce a primal direction for var getting all of its data
 void process_variable_dir( ColVariable & var );

 /// produce a primal direction to a  "linear" variable
 void process_var_dir_l( ColVariable & var , VarValue l , VarValue u ,
			 OFValue b );

 /// produce a primal direction to a  "quadratic" variable
 void process_var_dir_q( ColVariable & var , VarValue l , VarValue u ,
			 OFValue a );

 /// produce a dual scolution for var getting all of its data
 void process_variable_dual( ColVariable & var );

 /// produce a dual solution to a  "linear" variable
 void process_var_dual( ColVariable & var , VarValue l , VarValue u ,
			OFValue b ,
			OneVarConstraint * cl , OneVarConstraint * cu );

 /// produce a dual solution to a  "quadratic" variable
 void process_var_dual( ColVariable & var , VarValue l , VarValue u ,
			OFValue a , OFValue b ,
			OneVarConstraint * cl , OneVarConstraint * cu );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 char f_sol;        ///< whether the primal and/or dual solution is computed

 char f_sol_comp;   ///< the value of f_sol in the last call to compute()

 int f_state;       ///< the return status

 OFValue f_max_val;    ///< maximum of the Objective over the box

 OFValue f_min_val;    ///< minimum of the Objective over the box

 bool f_feas;       ///< if a feasible solution exists (and is loaded)

 char f_sense;
 ///< 1 if the Objective is max, 0 if it is min, -1 if it has to be computed

 Function * f_altobj;  ///< a possible alternative Function for the Objective

 Vec_Block f_desc;     ///< the Block and all its sub-Block, recursively

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };   // end( class BoxSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BoxSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File BoxSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
