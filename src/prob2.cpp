# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "prob2.h"

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r8_ceiling ( double r )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CEILING rounds an R8 "up" to the nearest integer.
//
//  Example:
//
//    R     R8_CEILING
//
//   -1.1  -1
//   -1.0  -1
//   -0.9   0
//    0.0   0
//    5.0   5
//    5.1   6
//    5.9   6
//    6.0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the real value to be rounded up.
//
//    Output, int R8_CEILING, the rounded value.
//
{
  int value;

  value = ( int ) ( r );

  if ( ( double ) value < r )
  {
    value = value + 1;
  }

  return value;
}
//****************************************************************************80

double r8_csc ( double theta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSC returns the cosecant of X.
//
//  Discussion:
//
//    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double THETA, the angle, in radians, whose cosecant is desired.
//    It must be the case that SIN ( THETA ) is not zero.
//
//    Output, double R8_CSC, the cosecant of THETA.
//
{
  double value;

  value = sin ( theta );

  if ( value == 0.0 )
  {
    //cout << " \n";
    //cout << "R8_CSC - Fatal error!\n";
    //cout << "  Cosecant undefined for THETA = " << theta << "\n";
    //exit ( 1 );
  }

  value = 1.0 / value;

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03,
    8.4171387781295E-04,
   -5.952379913043012E-04,
    7.93650793500350248E-04,
   -2.777777777777681622553E-03,
    8.333333333333333331554247E-02,
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01,
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02,
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04,
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02,
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03,
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03,
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" real value.
//
{
  return HUGE_VAL;
}
//****************************************************************************80


double bessel_i0 ( double arg )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0 evaluates the modified Bessel function I0.
//
//  Discussion:
//
//    The main computation evaluates slightly modified forms of
//    minimax approximations generated by Blair and Edwards, Chalk
//    River (Atomic Energy of Canada Limited) Report AECL-4928,
//    October, 1974.  This transportable program is patterned after
//    the machine dependent FUNPACK packet NATSI0, but cannot match
//    that version for efficiency or accuracy.  This version uses
//    rational functions that theoretically approximate I-SUB-0(X)
//    to at least 18 significant decimal digits.
//
//  Machine dependent constants:
//
//    beta   = Radix for the floating-point system
//    maxexp = Smallest power of beta that overflows
//    XSMALL = Positive argument such that 1.0 - X = 1.0 to
//             machine precision for all ABS(X) .LE. XSMALL.
//    XMAX =   Largest argument acceptable to BESI0;  Solution to
//             equation:
//               W(X) * (1+1/(8*X)+9/(128*X^2) = beta**maxexp
//             where  W(X) = EXP(X)/sqrt(2*PI*X)
//
//    Approximate values for some important machines are:
//
//                             beta       maxexp       XSMALL
//
//    CRAY-1        (S.P.)       2         8191       3.55D-15
//    Cyber 180/855
//      under NOS   (S.P.)       2         1070       3.55D-15
//    IEEE (IBM/XT,
//      SUN, etc.)  (S.P.)       2          128       2.98D-8
//    IEEE (IBM/XT,
//      SUN, etc.)  (D.P.)       2         1024       5.55D-17
//    IBM 3033      (D.P.)      16           63       6.95D-18
//    VAX           (S.P.)       2          127       2.98D-8
//    VAX D-Format  (D.P.)       2          127       6.95D-18
//    VAX G-Format  (D.P.)       2         1023       5.55D-17
//
//
//                                  XMAX
//
//    CRAY-1        (S.P.)       5682.810
//    Cyber 180/855
//      under NOS   (S.P.)       745.893
//    IEEE (IBM/XT,
//      SUN, etc.)  (S.P.)        91.900
//    IEEE (IBM/XT,
//      SUN, etc.)  (D.P.)       713.986
//    IBM 3033      (D.P.)       178.182
//    VAX           (S.P.)        91.203
//    VAX D-Format  (D.P.)        91.203
//    VAX G-Format  (D.P.)       713.293
//
//  Author:
//
//    Original FORTRAN77 version by W. J. Cody and L. Stoltz.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ARG, the argument.
//
//    Output, double BESSEL_I0, the value of the modified Bessel function
//    of the first kind.
//
{
  double a;
  double b;
  double exp40 = 2.353852668370199854E+17;
  int i;
  double p[15] = {
    -5.2487866627945699800E-18,
    -1.5982226675653184646E-14,
    -2.6843448573468483278E-11,
    -3.0517226450451067446E-08,
    -2.5172644670688975051E-05,
    -1.5453977791786851041E-02,
    -7.0935347449210549190,
    -2.4125195876041896775E+03,
    -5.9545626019847898221E+05,
    -1.0313066708737980747E+08,
    -1.1912746104985237192E+10,
    -8.4925101247114157499E+11,
    -3.2940087627407749166E+13,
    -5.5050369673018427753E+14,
    -2.2335582639474375249E+15 };
  double pp[8] = {
    -3.9843750000000000000E-01,
     2.9205384596336793945,
    -2.4708469169133954315,
     4.7914889422856814203E-01,
    -3.7384991926068969150E-03,
    -2.6801520353328635310E-03,
     9.9168777670983678974E-05,
    -2.1877128189032726730E-06 };
  double q[5] = {
    -3.7277560179962773046E+03,
     6.5158506418655165707E+06,
    -6.5626560740833869295E+09,
     3.7604188704092954661E+12,
    -9.7087946179594019126E+14 };
  double qq[7] = {
    -3.1446690275135491500E+01,
     8.5539563258012929600E+01,
    -6.0228002066743340583E+01,
     1.3982595353892851542E+01,
    -1.1151759188741312645,
     3.2547697594819615062E-02,
    -5.5194330231005480228E-04 };
  double rec15 = 6.6666666666666666666E-02;
  double sump;
  double sumq;
  double value;
  double x;
  double xmax = 91.9;
  double xsmall = 2.98E-08;
  double xx;

  x = r8_abs ( arg );

  if ( x < xsmall )
  {
    value = 1.0;
  }
  else if ( x < 15.0 )
  {
//
//  XSMALL <= ABS(ARG) < 15.0
//
    xx = x * x;
    sump = p[0];
    for ( i = 1; i < 15; i++ )
    {
      sump = sump * xx + p[i];
    }

    xx = xx - 225.0;
    sumq = ((((
           xx + q[0] )
         * xx + q[1] )
         * xx + q[2] )
         * xx + q[3] )
         * xx + q[4];

    value = sump / sumq;
  }
  else if ( 15.0 <= x )
  {
    if ( xmax < x )
    {
      value = r8_huge ( );
    }
    else
    {
//
//  15.0 <= ABS(ARG)
//
      xx = 1.0 / x - rec15;

      sump = ((((((
                  pp[0]
           * xx + pp[1] )
           * xx + pp[2] )
           * xx + pp[3] )
           * xx + pp[4] )
           * xx + pp[5] )
           * xx + pp[6] )
           * xx + pp[7];

      sumq = ((((((
             xx + qq[0] )
           * xx + qq[1] )
           * xx + qq[2] )
           * xx + qq[3] )
           * xx + qq[4] )
           * xx + qq[5] )
           * xx + qq[6];

      value = sump / sumq;
//
//  Calculation reformulated to avoid premature overflow.
//
      if ( x <= xmax - 15.0 )
      {
        a = exp ( x );
        b = 1.0;
      }
      else
      {
        a = exp ( x - 40.0 );
        b = exp40;
      }

      value = ( ( value * a - pp[0] * a ) / sqrt ( x ) ) * b;
    }
  }

  return value;
}
//****************************************************************************80

double von_mises_pdf ( double x, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_PDF evaluates the von Mises PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
//
//    where:
//
//      I0(*) is the modified Bessel function of the first
//      kind of order 0.
//
//    The von Mises distribution for points on the unit circle is
//    analogous to the normal distribution of points on a line.
//    The variable X is interpreted as a deviation from the angle A,
//    with B controlling the amount of dispersion.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 160.
//
//    D J Best, N I Fisher,
//    Efficient Simulation of the von Mises Distribution,
//    Applied Statistics,
//    Volume 28, Number 2, pages 152-157.
//
//    Kanti Mardia, Peter Jupp,
//    Directional Statistics,
//    Wiley, 2000, QA276.M335
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    A - PI <= X <= A + PI.
//
//    Input, double A, B, the parameters of the PDF.
//    -PI <= A <= PI,
//    0.0 < B.
//
//    Output, double VON_MISES_PDF, the value of the PDF.
//
{
  double pdf;
  const double pi = 3.14159265358979323;

  if ( x < a - pi )
  {
    pdf = 0.0;
  }
  else if ( x <= a + pi )
  {
    pdf = exp ( b * cos ( x - a ) ) / ( 2.0 * pi * bessel_i0 ( b ) );
  }
  else
  {
    pdf = 0.0;
  }

  return pdf;
}
//****************************************************************************80