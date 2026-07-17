/*--------------------------------------------------------------------------*/
/*------------------------- File DQuadFunction.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DQuadFunction, which implements
 * C15Function with a diagonal (separable) quadratic function.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DQuadFunction
#define __DQuadFunction
                     /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "C15Function.h"

#include "ColVariable.h"

#include "Observer.h"

#include <cmath>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DQuadFunction ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a diagonal quadratic Function
/**< The class DQuadFunction implements C15Function with a diagonal
 * quadratic function of the form
 * \f[
 * f(x) = c + sum{i = 1, ..., n} ( a_i * x_i^2 + b_i * x_i )
 * \f]
 * where the scalar c is the constant term of the Function, and a_i and b_i
 * are the coefficients of the Variable x_i in the quadratic and linear
 * terms, respectively.
 *
 * This Function issues the following :Modification:
 *
 * - A C05FunctionModVarsAddd when Variable are added; note that, being this
 *   a separable Function, the addition is strongly quasi-additive.
 *
 * - A C05FunctionModVarsRngd or C05FunctionModVarsSbst when Variable are
 *   removed (depending if they were in a range or not); note that, being
 *   this a separable Function, the removal is strongly quasi-additive.
 *
 * - a C05FunctionModLinRngd or a C05FunctionModLinSbst when the *linear*
 *   coefficients *only* of some Variable change (depending if they were in
 *   a range or not); the shift is NaN, no check on the sign of the
 *   Variables / coefficient change is performed (yet).
 *
 * - a C05FunctionModRngd or a C05FunctionModSbst when the *quadratic*
 *   coefficients (and, possibly the linear ones as well) of some Variable
 *   change, (depending if they were in a range or not); the shift is NaN,
 *   no check on the sign of the coefficients is performed (yet) [note that
 *   it is somewhat easier to prove increase/decrease only for a quadratic
 *   function, think e.g. of the case of no linear coefficients, all
 *   positive a_i which all grow].
 *
 * - a FunctionMod when the constant term changes, with the shift equal to
 *   the difference between the new and old constant term values.
 *
 * TODO: when "a lot" of the linear and/or quadratic coefficients change,
 *       rather issue a C05FunctionMod of type AllEntriesChanged, because it
 *       might be faster to just re-read all the coefficients than to skip
 *       the few non-changed ones. */

class DQuadFunction : public C15Function {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 /// type of the coefficients of linear and quadratic terms = FunctionValue
 using Coefficient = FunctionValue;

 using v_coeff = std::vector< Coefficient >;  ///< a vector of Coefficients

 using v_coeff_it = v_coeff::iterator;  ///< iterator in v_coeff

 using c_v_coeff_it = v_coeff::const_iterator;  ///< const iterator in v_coeff

 /// one term of the sum: < ColVariable * , Coefficient , Coefficient >
 /** Triple < ColVariable * , Coefficient , Coefficient > that describes
  * the term of the sum corresponding to one [Col]Variable. */
 using coeff_triple = std::tuple< ColVariable * , Coefficient , Coefficient >;

 using c_coeff_triple = const coeff_triple;  ///< a const coeff_triple
 
 /// a vector of coeff_triple
 using v_coeff_triple = std::vector< coeff_triple >;

 /// a const vector of coeff_triple
 using v_c_coeff_triple = const v_coeff_triple;

/*--------------------------------------------------------------------------*/
 /// virtualized concrete iterator
 /** A concrete class deriving from ThinVarDepInterface::v_iterator and
  * implementing the concrete iterator for sifting through the "active"
  * Variable of a DQuadFunction. */

 class v_iterator : public ThinVarDepInterface::v_iterator {
  public:

  explicit v_iterator( v_coeff_triple::iterator & itr ) : itr_( itr ) {}
  explicit v_iterator( v_coeff_triple::iterator && itr )
   : itr_( std::move( itr ) ) {}

  v_iterator * clone( void ) override final {
   return( new v_iterator( itr_ ) );
   }

  void operator++( void ) override final { ++(itr_); }

  reference operator*( void ) const override final {
   return( *std::get< 0 >( *itr_ ) );
   }

  pointer operator->( void ) const override final {
   return( std::get< 0 >( *itr_ ) );
   }

  bool operator==( const ThinVarDepInterface::v_iterator & rhs )
   const override final {
   #ifdef NDEBUG
    auto tmp = static_cast< const DQuadFunction::v_iterator * >( & rhs );
    return( itr_ == tmp->itr_ );
   #else
    auto tmp = dynamic_cast< const DQuadFunction::v_iterator * >( &rhs );
    return( tmp ? itr_ == tmp->itr_ : false );
   #endif
   }

  bool operator!=( const ThinVarDepInterface::v_iterator & rhs )
   const override final {
   #ifdef NDEBUG
    auto tmp = static_cast< const DQuadFunction::v_iterator * >( & rhs );
    return( itr_ != tmp->itr_ );
   #else
    auto tmp = dynamic_cast< const DQuadFunction::v_iterator * >( &rhs );
    return( tmp ? itr_ != tmp->itr_ : true );
   #endif
   }

  private:

  v_coeff_triple::iterator itr_;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// virtualized concrete const_iterator
 /** A concrete class deriving from ThinVarDepInterface::v_const_iterator and
  * implementing the concrete iterator for sifting through the "active"
  * Variable of a DQuadFunction. */

 class v_const_iterator : public ThinVarDepInterface::v_const_iterator {
  public:

  explicit v_const_iterator( v_c_coeff_triple::const_iterator & itr )
   : itr_( itr ) {}
  explicit v_const_iterator( v_c_coeff_triple::const_iterator && itr )
   : itr_( std::move( itr ) ) {}

  v_const_iterator * clone( void ) override final {
   return( new v_const_iterator( itr_ ) );
   }

  void operator++( void ) override final { ++(itr_); }

  reference operator*( void ) const override final {
   return( *std::get< 0 >( *itr_ ) );
   }

  pointer operator->( void ) const override final {
   return( std::get< 0 >( *itr_ ) );
   }

  bool operator==( const ThinVarDepInterface::v_const_iterator & rhs )
   const override final {
   #ifdef NDEBUG
    auto tmp = static_cast< const DQuadFunction::v_const_iterator * >( & rhs );
    return( itr_ == tmp->itr_ );
   #else
    auto tmp = dynamic_cast< const DQuadFunction::v_const_iterator * >( &rhs );
    return( tmp ? itr_ == tmp->itr_ : false );
   #endif
   }

  bool operator!=( const ThinVarDepInterface::v_const_iterator & rhs )
   const override final {
   #ifdef NDEBUG
    auto tmp = static_cast< const DQuadFunction::v_const_iterator * >( & rhs );
    return( itr_ != tmp->itr_ );
   #else
    auto tmp = dynamic_cast< const DQuadFunction::v_const_iterator * >( &rhs );
    return( tmp ? itr_ != tmp->itr_ : true );
   #endif
   }

  private:

  v_coeff_triple::const_iterator itr_;
  };

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of DQuadFunction, taking the data describing it
 /** Constructor of DQuadFunction. It accepts:
  *
  * @param v_var, a && to a vector of triples < ColVariable * , Coefficient ,
  *        Coefficient >, with the first coefficient being that of the linear
  *        term and the second being that of the quadratic term of the
  *        function corresponding to the given ColVariable; as the &&
  *        tells, vars is "consumed" by the constructor and its resources
  *        become property of the LinearFunction object.
  *
  * @param ct, a FunctionValue providing the value of the constant term of the
  *        function;
  *
  * @param observer, a pointer to the Observer of this DQuadFunction.
  *
  * All inputs have a default ({}, 0, and nullptr, respectively) so that this
  * can be used as the void constructor.
  *
  * Important note:
  *
  *     THE ORDER OF THE vars VECTOR WILL DICTATE THE ORDER OF THE "ACTIVE"
  *     [Col]Variable OF THE DQuadFunction
  *
  * That is, get_active_var( 0 ) == std::get< 0 >( vars[ 0 ] ),
  * get_active_var( 1 ) == std::get< 0 >( vars[ 1 ] ), ... */

 explicit DQuadFunction( v_coeff_triple && v_var = {} ,
                         const FunctionValue ct = 0 ,
                         Observer * const observer = nullptr )
  : C15Function( observer ) , v_triples( std::move( v_var ) ) ,
    f_value( Inf< FunctionValue >() ), f_constant_term( ct ) { }

/*--------------------------------------------------------------------------*/
 /// destructor: it does nothing (explicitly)

 ~DQuadFunction() override = default;

/*--------------------------------------------------------------------------*/
 /// clear method: clears the v_triples field
 /** Method to "clear" the DQuadFunction: it clear() the vector v_triples.
  * This destroys the list of "active" Variable without unregistering from
  * them. Not that the LinearFunction would have, but the Observer may.
  * By not having any Variable, the Observer can no longer do that. */

 void clear( void ) override { v_triples.clear(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the whole (empty) set of parameters in one blow
 /** Although a DQuadFunction formally has a lot of parameters, in fact it
  * "listens to no-one"; hence, the implementation of set_ComputeConfig() is
  * almost (but not quite) a trivial one. The issue is that DQuadFunction
  * does not currently implement the global pool management, and therefore
  * it will throw exception if a nonzero global pool size is set. */

 void set_ComputeConfig( ComputeConfig * scfg ) override final
 {
  if( ! scfg )  // setting nothing
   return;      // nothing to do

  for( auto pair : scfg->int_pars )
   if( int_par_str2idx( pair.first ) == intGPMaxSz ) {
    if( pair.second > 0 )
     throw( std::invalid_argument( 
			"global pool not supported yet by DQuadFunction" ) );
    break;
    }
  }

/*--------------------------------------------------------------------------*/
 /// set a given integer (int) numerical parameter
 /** Set a given integer (int) numerical parameter. DQuadFunction in fact
  * "listens to no parameter"; hence, the implementation of set_par() should
  * be void. However, DQuadFunction does not currently implement the global
  * pool management, and therefore it will throw exception if a nonzero
  * global pool size is set. */

 /* FIX (overload hiding): declaring any set_par() here hides ALL base-class
  * set_par() overloads. With no double overload and no using-declaration, a
  * set_par( par , <double> ) call narrows the double to int and dispatches to
  * the int overload - whose enum is a SEPARATE index space - silently
  * corrupting an unrelated parameter. Same defect as MILPSolver, where
  * dblAAccSol (double idx 7) overwrote intRelaxIntVars (int idx 7). */
 using C15Function::set_par;
 void set_par( const idx_type par , const int value ) override final
 {
  if( ( par == intGPMaxSz ) && ( value > 0 ) )
   throw( std::invalid_argument( 
			"global pool not supported yet by DQuadFunction" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE DQuadFunction ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the DQuadFunction
    @{ */

 /// returns the vector of triples (ColVariable *, Coefficient, Coefficient)

 v_c_coeff_triple & get_v_var( void ) const { return( v_triples ); }

/*--------------------------------------------------------------------------*/
 /// returns the linear Coefficient of the i-th Variable
 /** This method returns the linear Coefficient of the i-th Variable of this
  * DQuadFunction. The index i must be between 0 and get_num_active_var() - 1.
  *
  * @param i Index of the Variable whose linear coefficient is desired. */

 Coefficient get_linear_coefficient( Index i ) const {
  return( std::get< 1 >( *( v_triples.begin() + i ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the quadratic Coefficient of the i-th Variable
 /** This method returns the linear Coefficient of the i-th Variable of this
  * DQuadFunction. The index i must be between 0 and get_num_active_var() - 1.
  *
  * @param i Index of the Variable whose quadratic coefficient is desired. */

 Coefficient get_quadratic_coefficient( Index i ) const {
  return( std::get< 2 >( *( v_triples.begin() + i ) ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS DESCRIBING THE BEHAVIOR OF THE DQuadFunction -----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods describing the behavior of the DQuadFunction
 *  @{ */

 int compute( bool changedvars ) override;

/*--------------------------------------------------------------------------*/
 /// returns the value of the DQuadFunction

 FunctionValue get_value( void ) const override { return( f_value ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the DQuadFunction is exact, hence lower_estimate == value

 FunctionValue get_lower_estimate( void ) const override final {
  return( f_value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the DQuadFunction is exact, hence upper_estimate == value

 FunctionValue get_upper_estimate( void ) const override final {
  return( f_value );
  }

/*--------------------------------------------------------------------------*/

 bool is_convex( void ) const override;

/*--------------------------------------------------------------------------*/

 bool is_concave( void ) const override;

/*--------------------------------------------------------------------------*/

 bool is_lower_semicontinuous( void ) const override final { return( true ); }

/*--------------------------------------------------------------------------*/

 bool is_upper_semicontinuous( void ) const override final { return( true ); }

/*--------------------------------------------------------------------------*/

 bool is_linear( void ) const override;

/*--------------------------------------------------------------------------*/

 void compute_hessian_approximation( void ) override { };

/*--------------------------------------------------------------------------*/

 void get_hessian_approximation( DenseHessian & hessian ) const override;

/*--------------------------------------------------------------------------*/

 void get_hessian_approximation( SparseHessian & hessian ) const override;

/*--------------------------------------------------------------------------*/

 bool is_continuously_differentiable( void ) const override final {
  return( true );
  }

/*--------------------------------------------------------------------------*/

 bool is_twice_continuously_differentiable( void ) const override final {
  return( true );
  }

/*--------------------------------------------------------------------------*/
 /// returns the linearization coefficient of the i-th active Variable

 double get_linearization_coefficient( c_Index i ) {
  return( 2 * std::get<0>( v_triples[ i ] )->get_value() *
	  std::get<2>( v_triples[ i ] ) + std::get<1>( v_triples[ i ] ) );
  }

/*--------------------------------------------------------------------------*/

 void get_linearization_coefficients( FunctionValue * g ,
				      Range range = INFRange ,
				      Index name = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_linearization_coefficients( SparseVector & g ,
				      Range range = INFRange ,
				      Index name = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/

 void get_linearization_coefficients( FunctionValue * g , c_Subset & subset ,
				      const bool ordered = false ,
				      Index name = Inf< Index >() ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void get_linearization_coefficients( SparseVector & g , c_Subset & subset ,
				      const bool ordered = false ,
				      Index name = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/
 /// returns the linearization constant of the current point
 /** Since DQuadFunction currently does not support the global pool, the
  * method only works for returning the value of the linearization constant
  * at the current point x, which is given by c - x^T A x. */

 FunctionValue get_linearization_constant( Index name = Inf< Index >() )
  override final {
  if( name < Inf< Index >() )
   throw( std::invalid_argument(
			"global pool not supported yet by DQuadFunction" ) );

  FunctionValue quadratic_term = 0.0;

  for( const auto & triple : v_triples ) {
   auto variable_value = std::get< 0 >( triple )->get_value();
   auto quadratic_coefficient_value = std::get< 2 >( triple );
   quadratic_term += variable_value * variable_value *
                     quadratic_coefficient_value;
   }

 return( this->f_constant_term - quadratic_term );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 ///< Returns the value of the constant term of this DQuadFunction.

 [[nodiscard]] FunctionValue get_constant_term( void ) const override {
  return( f_constant_term );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the DQuadFunction
 *  @{ */

 ///< get the whole (empty) set of parameters in one blow
 /** Although a DQuadFunction formally has a lot of parameters, in fact it
  * "listens to no-one"; hence, the implementation of get_ComputeConfig() is
  * quite a trivial one. */

 ComputeConfig * get_ComputeConfig( bool all , ComputeConfig * ocfg )
  const override final { return( nullptr ); }

/** @} ---------------------------------------------------------------------*/
/*----- METHODS FOR HANDLING "ACTIVE" Variable IN THE DQuadFunction --------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling the set of "active" Variable in the
 * DQuadFunction; this is the actual concrete implementation exploiting the
 * vector v_triples of triples.
 * @{ */

 Index get_num_active_var() const override final { return( v_triples.size() ); }

/*--------------------------------------------------------------------------*/

 Index is_active( const Variable * const var ) const override final
 {
  auto idx = std::find_if( v_triples.begin() , v_triples.end() ,
			   [ & var ]( const auto & p ) -> bool {
			    return( std::get< 0 >( p ) == var );
			    } );
  return( idx != v_triples.end() ? std::distance( v_triples.begin(), idx )
	                         : Inf< Index >() );
  }

/*--------------------------------------------------------------------------*/

 void map_active( c_Vec_p_Var & vars , Subset & map , bool ordered )
  const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 Subset map_index( const std::vector< Variable * > & vars , c_Subset & nms )
  const override {
  if( vars.size() != nms.size() )
   throw( std::invalid_argument( "vars and nms sizes do not match" ) );

  Subset map( nms.size() );
  if( map.empty() )
   return( map );

  auto nmsit = nms.begin();
  auto varsit = vars.begin();
  auto mapit = map.begin();

  // for all Variable in the set
  while( varsit != vars.end() ) {
   auto var = *(varsit++);  // next variable
   auto oi = *(nmsit++);    // its original index
   // if var has not been deleted (and, possibly, re-added), its index
   // must be <= oi: search backward from oi to find it
   auto avoi = std::get< 0 >( v_triples[ oi ] );
   Index i = oi;
   while( var != avoi ) {
    if( ! i )
     break;
    avoi = std::get< 0 >( v_triples[ --i ] );
    }
   if( var == avoi )  // the Variable was found
    *(mapit++) = i;   // this is its index
   else {             // the Variable was not found
    // restart the search from the last variable to oi (excluded), for the
    // case where var has been deleted and re-added, and therefore its
    // index can now be arbitrary (but it is more likely to be "close to
    // the end" than "at the beginning")
    for( i = v_triples.size() , avoi = std::get< 0 >( v_triples[ --i ] ) ;
	 ( var != avoi ) && ( i > oi ) ; )
     avoi = std::get< 0 >( v_triples[ --i ] );
    *(mapit++) = ( var == avoi ) ? i : Inf< Index >();
    }
   }  // end( for all Variable )

  return( map );

  }  // end( map_index( Subset )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 Subset map_index( const std::vector< Variable * > & vars , Range rng )
  const override {
  if( vars.size() != rng.second - rng.first )
   throw( std::invalid_argument( "vars and rng sizes do not match" ) );

  Subset map( vars.size() );
  if( map.empty() )
   return( map );

  auto varsit = vars.begin();
  auto mapit = map.begin();

  // for all Variable in the set
  for( auto oi = rng.first ; oi < rng.second ; ++oi ) {
   auto var = *(varsit++);  // next variable
   // if var has not been deleted (and, possibly, re-added), its index
   // must be <= oi: search backward from oi to find it
   auto avoi = std::get< 0 >( v_triples[ oi ] );
   Index i = oi;
   while( var != avoi ) {
    if( ! i )
     break;
    avoi = std::get< 0 >( v_triples[ --i ] );
    }
   if( var == avoi )  // the Variable was found
    *(mapit++) = i;   // this is its index
   else {             // the Variable was not found
    // restart the search from the last variable to oi (excluded), for the
    // case where var has been deleted and re-added, and therefore its
    // index can now be arbitrary (but it is more likely to be "close to
    // the end" than "at the beginning")
    for( i = v_triples.size() , avoi = std::get< 0 >( v_triples[ --i ] ) ;
	 ( var != avoi ) && ( i > oi ) ; )
     avoi = std::get< 0 >( v_triples[ --i ] );
    *(mapit++) = ( var == avoi ) ? i : Inf< Index >();
    }
   }  // end( for all Variable )

  return( map );

  }  // end( map_index( Range )

/*--------------------------------------------------------------------------*/

 Variable * get_active_var( const Index i ) const override final {
  return( std::get< 0 >( v_triples[ i ] ) );
  }

/*--------------------------------------------------------------------------*/

 v_iterator * v_begin( void ) override final {
  return( new DQuadFunction::v_iterator( v_triples.begin() ) );
  }

/*--------------------------------------------------------------------------*/

 v_const_iterator * v_begin( void ) const override final {
  return( new DQuadFunction::v_const_iterator( v_triples.begin() ) );
  }

/*--------------------------------------------------------------------------*/

 v_iterator * v_end( void ) override final {
  return( new DQuadFunction::v_iterator( v_triples.end() ) );
  }

/*--------------------------------------------------------------------------*/

 v_const_iterator * v_end( void ) const override final {
  return( new DQuadFunction::v_const_iterator( v_triples.end() ) );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR MODIFYING THE DQuadFunction -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DQuadFunction
 *  @{ */

 /// add a set of new Variable to the DQuadFunction
 /**< Method that receives a pointer to a vector of triples < ColVariable * ,
  * Coefficient , Coefficient > and adds these triples to the list of
  * triples already in the DQuadFunction. The first coefficient is that in
  * the linear term and the second one is that in the quadratic term. Note
  * that
  *
  *     IT IS EXPECTED THAT NONE OF THE NEW ColVariable IS ALREADY "ACTIVE"
  *     IN THE DQuadFunction, BUT NO CHECK IS DONE TO ENSURE THIS
  *
  * Indeed, the check is costly, and the DQuadFunction does not really have
  * a functional issue with repeated ColVariable. The real issue rather comes
  * whenever the DQuadFunction is used within a Constraint or Objective that
  * need to register itself among the "active" Variable of the DQuadFunction;
  * this process is not structured to work with multiple copies of the same
  * "active" Variable. Thus, a DQuadFunction used within such an object
  * should not have repeated Variable, but if this is an issue then the
  * check will have to be performed elsewhere, it is not done here.
  *
  * As the && tells, vars is "consumed" by the method and its resources
  * become property of the DQuadFunction object.
  *
  * The parameter issueMod decides if and how the C05FunctionModVars is
  * issued, as described in Observer::make_par(). Note that a diagonal
  * quadratic function is additive, and therefore strongly quasi-additive. */

 void add_variables( v_coeff_triple && vars ,
                     ModParam issueMod = eModBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// add one single new Variable to the DQuadFunction
 /** Like add_variables(), but just only one Variable. The linear_coeff is
  * the coefficient of the Variable in the linear term and quadratic_coeff
  * is the coefficient of the Variable in the quadratic term. */

 void add_variable( ColVariable * var , Coefficient lin_coeff ,
                    Coefficient quad_coeff , ModParam issueMod = eModBlck );

/*--------------------------------------------------------------------------*/
 /// modify a single existing quadratic term
 /** Method that modifies both the linear and the quadratic coefficients for
  * the i-th "active" Variable; if i is not a valid index, exception is
  * thrown. 
  *
  * The parameter issueMod decides if and how the C05FunctionModRngd is
  * issued, as described in Observer::make_par(). */

 void modify_term( Index i , Coefficient lin_coeff ,
                   Coefficient quad_coeff , ModParam issueMod = eModBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// modify a single existing linear coefficient
 /** Method that modifies only the linear coefficient for the i-th "active"
  * Variable; if i is not a valid index, exception is thrown.
  *
  * The parameter issueMod decides if and how the C05FunctionModLinRngd is
  * issued, as described in Observer::make_par(). This is precisely the
  * advantage of  using this method w.r.t. modify_term( i , lin , quad )
  * with quad == 0; a "less general" C05FunctionModLinRngd can be issued in
  * place of a C05FunctionModRngd one. */

 void modify_linear_coefficient( Index i , Coefficient coeff ,
                                 ModParam issueMod = eModBlck );

/*--------------------------------------------------------------------------*/
 /// modify a set of existing quadratic terms
 /** Method that receives (an iterator into) the vector of quadratic
  * coefficient values, (an iterator into) the vector of linea coefficient
  * values, and the set of index of the ColVariable whose (both) coefficients
  * need be modified, and sets the quadratic coefficient of ColVariable
  * nms[ i ] to *( NQuadCoef + i ), and the linear one to *( NLinCoef + i ).
  * As the && tells, nms becomes property of the DQuadFunction object
  * (possibly to be immediately dispatched to the issued C05FunctionModSbst).
  * The ordered parameter tells if subset is ordered by increasing Index,
  * which is actually *not* helpful for DQuadFunction but it may be for 
  * Block/Solver having to deal with the C05FunctionModSbst.

  * The parameter issueMod decides if and how the C05FunctionModSbst is
  * issued, as described in Observer::make_par(). */

 void modify_terms( c_v_coeff_it NQuadCoef , c_v_coeff_it NLinCoef ,
                    Subset && nms , bool ordered = false ,
                    ModParam issueMod = eModBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// modify a set of existing linear coefficients
 /** Method that receives the vector of linear coefficient values and the set
  * of indices of the ColVariable whose linear coefficients need be modified,
  * and sets the linear coefficient of ColVariable nms[ i ] to NCoef[ i ]. As
  * the && tells, both NCoef and nms become property of the DQuadFunction
  * object (possibly to be immediately dispatched to the issued
  * C05FunctionModLinSbst). The ordered parameter tells if subset is ordered
  * by increasing Index, which is actually *not* helpful for DQuadFunction
  * but it may be for a Block/Solver having to deal with the
  * C05FunctionModLinSbst.
  *
  * The parameter issueMod decides if and how the C05FunctionModLinSbst is
  * issued, as described in Observer::make_par(). This is precisely the
  * advantage of using this method w.r.t. modify_terms() with all-zero
  * quadratic coefficients; a "less general" C05FunctionModLinSbst can be
  * issued in place of a C05FunctionModSbst one. */

 void modify_linear_coefficients( Vec_FunctionValue && NCoef ,
                                  Subset && nms , bool ordered = false ,
                                  ModParam issueMod = eModBlck );

/*--------------------------------------------------------------------------*/
 /// modify a range of existing quadratic terms
 /** Method that receives (an iterator into) the vector of quadratic
  * coefficient values, (an iterator into) the vector of linea coefficient
  * values, and modifies both the linear and quadratic coefficients of all
  * the Variable that are in position i, where range.first <= i <
  * min( range.second , get_num_active_var() ), giving them value 
  * respectively *( NQuadCoef + i - range.first ), and 
  * *( NLinCoef + i - range.first ).
  *
  * The parameter issueMod decides if and how the C05FunctionModSbst is
  * issued, as described in Observer::make_par(). */

 void modify_terms( c_v_coeff_it NQuadCoef , c_v_coeff_it NLinCoef ,
                    Range range = INFRange , ModParam issueMod = eModBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// modify a range of linear coefficients
 /** Modify the linear coefficients of all the Variable that are in position
  * i, where range.first <= i < min( range.second , get_num_active_var() ),
  * giving them value NCoef[ i - range.first ]. As the && tells, NCoef
  * becomes property of the LinearFunction object (possibly to be
  * immediately dispatched to the issued C05FunctionModLinRngd).
  *
  * The parameter issueMod decides if and how the C05FunctionModLinRngd is
  * issued, as described in Observer::make_par(). This is precisely the
  * advantage of using this method w.r.t. modify_terms() with all-zero
  * quadratic coefficients; a "less general" C05FunctionModLinRngd can be
  * issued in place of a C05FunctionModRngd one. */

 void modify_linear_coefficients( Vec_FunctionValue && NCoef ,
                                  Range range = INFRange ,
                                  ModParam issueMod = eModBlck );

/*--------------------------------------------------------------------------*/
 /// remove the i-th "active" Variable from the DQuadFunction
 /** Remove the i-th "active" Variable from the DQuadFunction. This is
  * *mathematically* equivalent to setting both its linear and quadratic
  * coefficients to zero, but it is considered a "stronger" operation (it is
  * possible to have an active Variable with both zero coefficients). If there
  * is no Variable with the given index, an exception is thrown.
  *
  * The parameter issueMod decides if and how the C05FunctionModVarsRngd is
  * issued, as described in Observer::make_par(). Note that a diagonal
  * quadratic function is additive, and therefore strongly quasi-additive,
  * which is why a C05FunctionModVarsRngd is issued as opposed to a
  * FunctionModVarsRngd one. */

 void remove_variable( Index i , ModParam issueMod = eModBlck )
  override final;

/*--------------------------------------------------------------------------*/
 /// remove a range of "active" Variable
 /** Remove all the "active" Variable in the given Range, i.e., all those
  * with index i s.t. range.first <= i < min( range.second ,
  * get_num_active_var() ), from this DQuadFunction.
  *
  * The parameter issueMod decides if and how the C05FunctionModVarsRngd is
  * issued, as described in Observer::make_par(). Note that a linear function
  * is additive, and therefore strongly quasi-additive, which is why a
  * C05FunctionModVarsRngd is issued as opposed to a FunctionModVarsRngd
  * one. */

 void remove_variables( Range range , ModParam issueMod = eModBlck )
  override final;
  
/*--------------------------------------------------------------------------*/
 /// remove the given subset of Variable
/** Remove all the Variable in the given set of indices. As the && tells,
  * nms becomes property of the DQuadFunction object (possibly to be
  * immediately dispatched to the issued C05FunctionModVarSbst). a special
  * setting is if
  *
  *     nms.empty() == true, IN WHICH CASE ALL Variable ARE ELIMINATED
  *
  * The parameter ordered tells if nms is ordered by increasing index. This
  * is useful for efficiently deleting them; indeed, if ordered == false the
  * vector is sorted inside. The parameter is provided to signal the lucky
  * case in which this operation can be avoided since nms "naturally" comes
  * out ordered.
  *
  * The parameter issueMod decides if and how the C05FunctionModVarSbst is
  * issued, as described in Observer::make_par(). Note that a linear function
  * is additive, and therefore strongly quasi-additive, which is why a
  * C05FunctionModVarSbst is issued as opposed to a FunctionModVarSbst one. */

 void remove_variables( Subset && nms , bool ordered = false ,
                        ModParam issueMod = eModBlck ) override final;

/*--------------------------------------------------------------------------*/
 ///< sets the value of the constant term of this function.
 /** Method that sets the new value to the constant term of this
  * diagonal quadratic Function to constant_term.
  *
  * The parameter issueMod decides if and how the C05FunctionMod (with
  * type() == NothingChanged and therefore clearly which().empty, and with
  * shift() == constant_term - < old constant term >) is issued, as
  * described in Observer::make_par(). */

 void set_constant_term( FunctionValue constant_term ,
                         ModParam issueMod = eModBlck );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// printing the DQuadFunction
 void print( std::ostream & output ) const override {
  output << "DQuadFunction [" << this << "] observed by ["
         << &f_Observer << "] with " << get_num_active_var()
         << " active variables;";
  output << "current value = " << get_value();
  }

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 v_coeff_triple v_triples;
 /**< vector of triples < ColVariable * , Coefficient , Coefficient >
  * characterizing the linear and the quadratic terms of the DQuadFunction;
  * the first coefficient is the linear one and the second one is the
  * quadratic term. */

 FunctionValue f_value;  ///< the value of the function

 FunctionValue f_constant_term = 0;
 ///< the value of the constant term of this function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

 };  // end( class( DQuadFunction ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* DQuadFunction.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File DQuadFunction.h ---------------------------*/
/*--------------------------------------------------------------------------*/
