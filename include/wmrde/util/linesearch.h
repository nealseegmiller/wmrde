#ifndef _WMRDE_LINESEARCH_H_
#define _WMRDE_LINESEARCH_H_

#include <wmrde/algebra/dynamic_matrix.h> //for Vecd

namespace wmrde
{

struct LinesearchOptions
{
  Real alpha_tol;
  Real c1;
  Real c2;

  LinesearchOptions() :
    alpha_tol(1e-3),
    c1(1e-4),
    c2(0.9)
  {}
};

/*!
 * perform line search. solve for step size (alpha) such that strong
 * Wolfe conditions are satisfied.
 * Reference:
 * Nocedal and Wright, Numerical Optimization, Springer, 1999. Algorithm 3.2 and 3.3
 * \param[in] p the search direction
 * \param[in] alpha_max the maximum step size
 * \param[in] fun function to compute cost and gradient. The signature must be:
 *            void fun(const Vecd& x, Real& cost, Vecd& gradient)
 * \param[in,out] x on return changed to x + alpha*p
 * \param[in,out] cost input cost at x, on return changed to value at x + alpha*p
 * \param[in,out] gradient input gradient at x, on return changed to value at x + alpha*p
 * \param[in] options tolerances for linesearch termination
 * \return the step size (alpha)
 */
template<typename FunctionType>
Real linesearch(const Vecd& p, const Real alpha_max, FunctionType fun,
    Vecd& x, Real& cost, Vecd& gradient,
    const LinesearchOptions options = LinesearchOptions())
{
  //backup values for alpha=0
  Real alpha=0; //for return
  Vecd x0 = x;
//  fun(x, cost, gradient); //assume already computed
  Real cost_0 = cost;
  Real grada_0 = gradient.dot(p); //The gradient of cost wrt alpha

  //first try taking the maximum step
  Real alpha_j = alpha_max;
  x += alpha_j * p;

  fun(x,cost,gradient);
  Real cost_j = cost;
  Real grada_j = gradient*p;

  bool dozoom = false;

  Real alpha_hi, alpha_lo, cost_lo;
  if (cost_j > cost_0 + options.c1*alpha_j*grada_0) {
    //decrease condition not satisfied
    dozoom = true;
    alpha_hi = alpha_j;
    alpha_lo = 0;
    cost_lo = cost_0;
  } else {
    if (std::abs(grada_j) <= -options.c2*grada_0) {
      //strong Wolfe conditions satisfied!
      alpha = alpha_j;
    } else {
      //curvature condition not satisfied
      if (grada_j >= 0) {
        dozoom = true;
        alpha_hi = 0;
        alpha_lo = alpha_j;
        cost_lo = cost_j;
      } else {
        //step could be larger
        alpha = alpha_j;
      }
    }
  }

  if (dozoom) {
    while (true) {
      if (fabs(alpha_hi-alpha_lo) < options.alpha_tol) {
        //close enough, break with alpha_lo
        if (alpha_j != alpha_lo) {
          x += alpha_lo*p;
          fun(x, cost, gradient); //update cost & gradient values
        }
        alpha = alpha_lo;
        break;
      }

      alpha_j = (alpha_lo + alpha_hi)/2; //bisection
      x += alpha_j*p;
      fun(x,cost,gradient);
      cost_j = cost;
      grada_j = gradient.dot(p);

      if (cost_j > cost_0 + options.c1*alpha_j*grada_0 || cost_j >= cost_lo) {
        //decrease condition not satisfied
        alpha_hi = alpha_j;
      } else {
        if (std::abs(grada_j) <= -options.c2*grada_0) {
          //strong Wolfe conditions satisfied!
          alpha = alpha_j;
          break;
        }
        if (grada_j*(alpha_hi-alpha_lo) >= 0) {
          alpha_hi = alpha_lo;
        }
        alpha_lo = alpha_j;
        cost_lo = cost_j;
      }
    }
  }

  return alpha;
}

} //namespace

/*

//line search, solve for step size (alpha) such that strong Wolfe conditions are satisfied
//n:			size of x
//p:			search direction, x =  x + alpha*p
//alpha_max:	maximum step size
//fCost:		cost function
//fGradient:	gradient function
//input initial values, output updated values for:
//x:		
//cost:
//grad:		gradient
//Reference:
//Nocedal and Wright, Numerical Optimization, Springer, 1999. Algorithm 3.2 and 3.3
template<typename Func1, typename Func2>
Real linesearch( const int n, const Real p[], const Real alpha_max, Func1 fCost, Func2 fGradient, 
				Real x[], Real& cost, Real grad[], ...) {

	const int MAXN = 20;
	Real alpha=0; //for return

	//options
	const Real alpha_tol = 1e-3;
	//for strong Wolfe conditions
	const Real c1 = 1e-4;
	const Real c2 = 0.9;

	
	//back up at alpha=0
	//grada: gradient of cost wrt alpha
	Real x0[MAXN]; //x (alpha = 0);
	Real cost_0 = cost;
	Real grada_0 = 0;
	for (int i=0; i<n; i++) {
		x0[i] = x[i]; //copy
		grada_0 += grad[i]*p[i]; //dot product
	}


	Real alpha_j, cost_j, grada_j;

	//first try taking the maximum step
	alpha_j = alpha_max;
	for (int i=0; i<n; i++) 
		x[i] = x[i] + alpha_j*p[i];

	cost = fCost(x);
	fGradient(grad);

	cost_j = cost;
	grada_j = 0; 
	for (int i=0; i<n; i++) 
		grada_j += grad[i]*p[i]; //dot product

	bool dozoom = false;

	Real alpha_hi, alpha_lo, cost_lo;
	if (cost_j > cost_0 + c1*alpha_j*grada_0) { 
		//decrease condition not satisfied
		dozoom = true;
		alpha_hi = alpha_j;
		alpha_lo = 0;
		cost_lo = cost_0;
	} else {
		//if (fabs(grada_j) <= c2*fabs(grada_0)) { //?
		if (fabs(grada_j) <= -c2*grada_0) {
			//strong Wolfe conditions satisfied!
			alpha = alpha_j;
		} else {
			//curvature condition not satisfied
			if (grada_j >= 0) {
				dozoom = true;
				alpha_hi = 0;
				alpha_lo = alpha_j;
				cost_lo = cost_j;
			} else {
				//step could be larger
				alpha = alpha_j;
			}
		}

	}

	if (dozoom) {
		while (true) {
			if (fabs(alpha_hi-alpha_lo) < alpha_tol) {
				//close enough, break with alpha_lo
				if (alpha_j != alpha_lo) {
					//need to call cost and gradient functions again
					//to set cost, grad, and va_list vars
					for (int i=0; i<n; i++) 
						x[i] = x0[i] + alpha_lo*p[i];
					cost = fCost(x);
					fGradient(grad);
				}
				alpha = alpha_lo;
				break;
			}

			alpha_j = (alpha_lo + alpha_hi)/2; //bisection
			for (int i=0; i<n; i++) 
				x[i] = x0[i] + alpha_j*p[i];
			cost = fCost(x);
			cost_j = cost;

			if (cost_j > cost_0 + c1*alpha_j*grada_0 || cost_j >= cost_lo) {
				//decrease condition not satisfied
				alpha_hi = alpha_j;
			} else {
				fGradient(grad);
				grada_j=0; 
				for (int i=0; i<n; i++) 
					grada_j += grad[i]*p[i]; //dot product
				if (fabs(grada_j) <= -c2*grada_0) {
					//strong Wolfe conditions satisfied!
					alpha = alpha_j;
					break;
				}
				if (grada_j*(alpha_hi-alpha_lo) >= 0) {
					alpha_hi = alpha_lo;
				}
				alpha_lo = alpha_j;
				cost_lo = cost_j;
			}
		}
	}

	return alpha;

}

*/

#endif //_WMRDE_LINESEARCH_H_
