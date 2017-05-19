#ifndef _WMRDE_ROOTFINDING_H_
#define _WMRDE_ROOTFINDING_H_

#include <wmrde/common.h>
#include <wmrde/rosout.h>

namespace wmrde
{

/*!
 * Find root via bisection.
 * http://en.wikipedia.org/wiki/Bisection_method
 * \param a lower limit for x value
 * \param b upper limit for x value
 * \param fx the function f(x)
 * \param tolx tolerance for x value
 * \param tolfx tolerance for f(x) value
 * \return x value for which f(x) is near zero. Return nan if fails
 */
template<typename Func>
Real findRootBisection(Real a, Real b, Func fx, const Real tolx, const Real tolfx)
{
  Real fa = fx(a);
  Real fb = fx(b);
  Real c,fc;

  if (fa*fb >= 0) //no root or multiple roots
  {
    return RealNan();
  }

  int iter = 0;
  while (iter < 20)
  {
    c = (a+b)/2;
    fc = fx(c);

    ROS_INFO("  findRootBisection() iter = %d: x = %f, f(x) = %f\n", iter, c, fc);

    if (fabs(fc) <= tolfx ||
        fabs(b-a)/2 <= tolx)
    {
      break; //within tolerance
    }

    if (fa*fc > 0) { a=c; fa=fc; }
    else { b=c; fb=fc; }

    iter++;
  }
  return c;
}

/*!
 * find root via Brent's method
 * http://en.wikipedia.org/wiki/Brent's_method
 * \param a lower limit for x value
 * \param b upper limit for x value
 * \param fx the function f(x)
 * \param tolx tolerance for x value
 * \param tolfx tolerance for f(x) value
 * \return x value for which f(x) is near zero. Return nan if fails
 */
template<typename Func>
Real findRootBrents( Real a, Real b, Func fx, const Real tolx, const Real tolfx)
{
  Real fa = fx(a);
  Real fb = fx(b);

  if (fa*fb >= 0) //no root or multiple roots
  {
    ROS_ERROR("No root or multiple roots");
    return RealNan();
  }

  if (fabs(fa) < fabs(fb))
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }

  Real c = a;
  Real fc = fa;
  Real d = RealInf();
  Real s;
  Real fs;

  bool mflag = true;
  int iter = 0;
  while ( iter < 20 )
  {
    if (fa != fc && fb != fc) {
      //inverse quadratic interpolation
      s = a*fb*fc/(fa-fb)/(fa-fc) + b*fa*fc/(fb-fa)/(fb-fc) + c*fa*fb/(fc-fa)/(fc-fb);
    } else {
      //secant rule
      s = b - fb*(b-a)/(fb-fa);
    }

    Real tmp = 3*a+b/4;
    bool cond1 = (s>tmp && s<b) || (s<tmp && s>b); cond1 = !cond1;
    bool cond2 =  mflag && fabs(s-b) >= fabs(b-c)/2;
    bool cond3 = !mflag && fabs(s-b) >= fabs(c-d)/2;
    bool cond4 =  mflag && fabs(b-c) < 2*tolx;
    bool cond5 = !mflag && fabs(c-d) < 2*tolx;
    //TODO, can cond4 or cond5 ever be true?

    if ( cond1 || cond2 || cond3 || cond4 || cond5 )
    {
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

    if (fabs(fa) < fabs(fb))
    {
      std::swap(a,b);
      std::swap(fa,fb);
    }

    ROS_INFO("  findRootBrents() iter = %d: x = %f, f(x) = %f\n", iter, b, fb);

    if (fabs(fb) <= tolfx ||
        fabs(a-b) <= 2*tolx)
    {
      break;
    }

    iter++;
  }

  return b;
}

} //namespace


#endif
