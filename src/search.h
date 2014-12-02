//search.h
//functions for root finding and line search
//use templates to pass in functions parameters
//alternative is to use std::function but that is much slower!

#ifndef _WMRSIM_SEARCH_H_
#define _WMRSIM_SEARCH_H_

//#include <functional> //for std::function, too slow!
#include <common/common.h>

//find root via bisection
//a:	lower bound
//b:	upper bound
//fx:	function, f(x)
//x:		tolerance on x
//tolfx:	tolerance on f(x)
template<typename Func>
Real findRootBisection( Real a, Real b, Func fx, const Real tolx, const Real tolfx) {
	//http://en.wikipedia.org/wiki/Bisection_method

	Real fa = fx(a);
	Real fb = fx(b);

	Real c,fc;

	if (fa*fb >= 0) {
		//TODO, no root or multiple roots?
		return REALNAN;
	}

	int iter = 0;

	while ( iter < 20 ) {

		c = (a+b)/2;
		fc = fx(c);

		//std::cout << "iter= " << iter << ", x= " << c << ", fx= " << fc << std::endl; //DEBUGGING

		if (fabs(fc) <= tolfx || fabs(b-a)/2 <= tolx) 
			break;

		if (fa*fc > 0) { a=c; fa=fc; } 
		else { b=c; fb=fc; }

		iter++;
		

	}
	return c;
}

//find root via Brent's method
//http://en.wikipedia.org/wiki/Brent's_method
template<typename Func>
Real findRootBrents( Real a, Real b, Func fx, const Real tolx, const Real tolfx) {
	

	Real fa = fx(a);
	Real fb = fx(b);

	if (fa*fb >= 0) {
		//TODO, no root or multiple roots?
		return REALNAN;
	}

	Real tmp;
	if (fabs(fa) < fabs(fb)) {
		//swap a,b
		tmp = a; a = b; b = tmp;
		tmp = fa; fa = fb; fb = tmp;
	}

	Real c = a;
	Real fc = fa;

	Real d = REALMAX;

	Real s;
	Real fs;

	Real tolx_ = 2*tolx;
	
	bool mflag = true;
	int iter = 0;
	

	bool cond1,cond2,cond3,cond4,cond5;

	while ( iter < 20 ) {

		
		if (fa != fc && fb != fc) {
			//inverse quadratic interpolation
			s = a*fb*fc/(fa-fb)/(fa-fc) + b*fa*fc/(fb-fa)/(fb-fc) + c*fa*fb/(fc-fa)/(fc-fb);
		} else {
			//secant rule
			s = b - fb*(b-a)/(fb-fa);
		}
		
		tmp = 3*a+b/4;
		cond1 = (s>tmp && s<b) || (s<tmp && s>b); cond1 = !cond1;
		cond2 =  mflag && fabs(s-b) >= fabs(b-c)/2;
		cond3 = !mflag && fabs(s-b) >= fabs(c-d)/2;
		cond4 =  mflag && fabs(b-c) < tolx_;
		cond5 = !mflag && fabs(c-d) < tolx_;
		//TODO, can cond4 or cond5 ever be true?

		if ( cond1 || cond2 || cond3 || cond4 || cond5 ) {
			//bisection method
			s=(a+b)/2;
			mflag=true;
		} else {
			mflag=false;
		}
		fs = fx(s);
		d = c;
		c = b;
		fc = fb;

		if (fa*fs < 0) { b=s; fb=fs; }
		else { a=s; fa=fs; }

		if (fabs(fa) < fabs(fb)) {
			//swap a,b
			tmp = a; a = b; b = tmp;
			tmp = fa; fa = fb; fb = tmp;
		}
		
		//std::cout << "iter= " << iter << ", x= " << b << ", fx= " << fb << std::endl; //DEBUGGING

		if ( fabs(fb) <= tolfx || fabs(a-b) <= tolx_ ) 
			break;

		iter++;
	}

	return b;
}


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


#endif 