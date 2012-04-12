//------------------------------------------------------------------------------
//     Main program of test problems in Report DAMTP 1992/NA5.
//------------------------------------------------------------------------------

using System;
using System.Linq;
using Cureos.Numerics;
using NUnit.Framework;

namespace Cscobyla.Tests
{
    [TestFixture]
    public class Cobyla2Tests
    {
        #region FIELDS

        private const int iprint = 3;
        private const double rhobeg = 0.5;

        #endregion

        #region TEST CASES

        [Test, Sequential]
        public void TestProblem1([Values(1.0e-6, 1.0e-8)] double rhoend, [Values(1.0e-5, 1.0e-7)] double acceptedError)
        {
            //     Minimization of a simple quadratic function of two variables.
            InvokeTestProblem(Calcfc1, 2, 0, rhoend, new[] { -1.0, 0.0 }, acceptedError);
        }

        public static void Calcfc1(int n, int m, double[] x, out double f, double[] con)
        {
            f = 10.0 * Math.Pow(x[0] + 1.0, 2.0) + Math.Pow(x[1], 2.0);
        }

        [Test, Sequential]
        public void TestProblem2([Values(1.0e-6, 1.0e-8)] double rhoend, [Values(1.0e-5, 1.0e-7)] double acceptedError)
        {
            //     Minimization of a simple quadratic function of two variables.
            InvokeTestProblem(Calcfc2, 2, 1, rhoend, new[] { Math.Sqrt(0.5), -Math.Sqrt(0.5) }, acceptedError);
        }

        public static void Calcfc2(int n, int m, double[] x, out double f, double[] con)
        {
            f = x[0] * x[1];
            con[0] = 1.0 - x[0] * x[0] - x[1] * x[1];
        }

        public void InvokeTestProblem(CalcfcDelegate calcfc, int n, int m, double rhoend, double[] xopt, double acceptedError)
        {
            var x = Enumerable.Repeat(1.0, n).ToArray();

            var maxfun = 3500;
            Cobyla2.Minimize(calcfc, n, m, x, rhobeg, rhoend, iprint, ref maxfun);

            var error = xopt.Zip(x, (xa, xb) => (xa - xb) * (xa - xb)).Sum();
            Assert.Less(error, acceptedError);
        }

        #endregion
    }
}
/*PROGRAM test_cobyla
USE common_nprob
USE cobyla2
IMPLICIT NONE

INTEGER, PARAMETER :: nn = 10

REAL (dp) :: rhobeg, rhoend, temp, tempa, tempb, tempc, tempd, x(nn),  &
             xopt(nn)
INTEGER   :: i, icase, iprint, m, maxfun, n

INTERFACE
  SUBROUTINE calcfc (n, m, x, f, con)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)
    INTEGER, INTENT(IN)    :: n, m
    REAL (dp), INTENT(IN)  :: x(:)
    REAL (dp), INTENT(OUT) :: f
    REAL (dp), INTENT(OUT) :: con(:)
  END SUBROUTINE calcfc
END INTERFACE

DO nprob=1,10
  IF (nprob == 1) THEN

//     Minimization of a simple quadratic function of two variables.

    WRITE(*, 10)
    10 FORMAT (/'       Output from test problem 1 (Simple quadratic)')
    n = 2
    m = 0
    xopt(1) = -1.0_dp
    xopt(2) = 0.0_dp
  ELSE IF (nprob == 2) THEN

//     Easy two dimensional minimization in unit circle.

    WRITE(*, 20)
    20 FORMAT (/'       Output from test problem 2 (2D unit circle calculation)')
    n = 2
    m = 1
    xopt(1) = SQRT(0.5_dp)
    xopt(2) = -xopt(1)
  ELSE IF (nprob == 3) THEN

//     Easy three dimensional minimization in ellipsoid.

    WRITE(*, 30)
    30 FORMAT (/'       Output from test problem 3 (3D ellipsoid calculation)')
    n = 3
    m = 1
    xopt(1) = 1.0_dp/SQRT(3.0_dp)
    xopt(2) = 1.0_dp/SQRT(6.0_dp)
    xopt(3) = -1.0_dp/3.0_dp
  ELSE IF (nprob == 4) THEN

//     Weak version of Rosenbrock's problem.

    WRITE(*, 40)
    40 FORMAT (/'       Output from test problem 4 (Weak Rosenbrock)')
    n = 2
    m = 0
    xopt(1) = -1.0_dp
    xopt(2) = 1.0_dp
  ELSE IF (nprob == 5) THEN

//     Intermediate version of Rosenbrock's problem.

    WRITE(*, 50)
    50 FORMAT (/'       Output from test problem 5 (Intermediate Rosenbrock)')
    n = 2
    m = 0
    xopt(1) = -1.0_dp
    xopt(2) = 1.0_dp
  ELSE IF (nprob == 6) THEN

//     This problem is taken from Fletcher's book Practical Methods of
//     Optimization and has the equation number (9.1.15).

    WRITE(*, 60)
    60 FORMAT (/'       Output from test problem 6 (Equation ',  &
               '(9.1.15) in Fletcher)')
    n = 2
    m = 2
    xopt(1) = SQRT(0.5_dp)
    xopt(2) = xopt(1)
  ELSE IF (nprob == 7) THEN

//     This problem is taken from Fletcher's book Practical Methods of
//     Optimization and has the equation number (14.4.2).

    WRITE(*, 70)
    70 FORMAT (/'       Output from test problem 7 (Equation ',  &
               '(14.4.2) in Fletcher)')
    n = 3
    m = 3
    xopt(1) = 0.0_dp
    xopt(2) = -3.0_dp
    xopt(3) = -3.0_dp
  ELSE IF (nprob == 8) THEN

//     This problem is taken from page 66 of Hock and Schittkowski's book Test
//     Examples for Nonlinear Programming Codes. It is their test problem Number
//     43, and has the name Rosen-Suzuki.

    WRITE(*, 80)
    80 FORMAT (/'       Output from test problem 8 (Rosen-Suzuki)')
    n=4
    m=3
    xopt(1) = 0.0_dp
    xopt(2) = 1.0_dp
    xopt(3) = 2.0_dp
    xopt(4) = -1.0_dp
  ELSE IF (nprob == 9) THEN

//     This problem is taken from page 111 of Hock and Schittkowski's
//     book Test Examples for Nonlinear Programming Codes. It is their
//     test problem Number 100.

    WRITE(*, 90)
    90 FORMAT (/'       Output from test problem 9 (Hock and Schittkowski 100)')
    n = 7
    m = 4
    xopt(1) = 2.330499_dp
    xopt(2) = 1.951372_dp
    xopt(3) = -0.4775414_dp
    xopt(4) = 4.365726_dp
    xopt(5) = -0.624487_dp
    xopt(6) = 1.038131_dp
    xopt(7) = 1.594227_dp
  ELSE IF (nprob == 10) THEN

//     This problem is taken from page 415 of Luenberger's book Applied
//     Nonlinear Programming. It is to maximize the area of a hexagon of
//     unit diameter.

    WRITE(*, 100)
    100 FORMAT (/'       Output from test problem 10 (Hexagon area)')
    n = 9
    m = 14
  END IF

  DO icase = 1,2
    x(1:n) = 1.0_dp
    rhobeg = 0.5_dp
    rhoend = 1.d-6
    IF (icase == 2) rhoend = 1.d-8
    iprint = 1
    maxfun = 3500
    CALL cobyla (n, m, x, rhobeg, rhoend, iprint, maxfun)
    IF (nprob == 10) THEN
      tempa = x(1) + x(3) + x(5) + x(7)
      tempb = x(2) + x(4) + x(6) + x(8)
      tempc = 0.5_dp/SQRT(tempa*tempa + tempb*tempb)
      tempd = tempc*SQRT(3.0_dp)
      xopt(1) = tempd*tempa + tempc*tempb
      xopt(2) = tempd*tempb - tempc*tempa
      xopt(3) = tempd*tempa - tempc*tempb
      xopt(4) = tempd*tempb + tempc*tempa
      DO i=1,4
        xopt(i+4) = xopt(i)
      END DO
    END IF
    temp = 0.0_dp
    DO i=1,n
      temp = temp + (x(i) - xopt(i))**2
    END DO
    WRITE(*, 150) SQRT(temp)
    150 FORMAT (/'     Least squares error in variables = ', G16.6)
  END DO
  WRITE(*, 170)
  170 FORMAT ('  ----------------------------------------------',  &
              '--------------------')
END DO

STOP
END PROGRAM test_cobyla


//------------------------------------------------------------------------------

SUBROUTINE calcfc (n, m, x, f, con)
USE common_nprob
IMPLICIT NONE

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: m
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: f
REAL (dp), INTENT(OUT)  :: con(:)

IF (nprob == 1) THEN

//     Test problem 1 (Simple quadratic)

  f = 10.0_dp*(x(1) + 1.0_dp)**2 + x(n)**2
ELSE IF (nprob == 2) THEN

//    Test problem 2 (2D unit circle calculation)

  f = x(1)*x(2)
  con(m) = 1.0_dp - x(1)**2 - x(n)**2
ELSE IF (nprob == 3) THEN

//     Test problem 3 (3D ellipsoid calculation)

  f = x(1)*x(2)*x(3)
  con(1) = 1.0_dp - x(1)**2 - 2.0_dp*x(2)**2 - 3.0_dp*x(n)**2
ELSE IF (nprob == 4) THEN

//     Test problem 4 (Weak Rosenbrock)

  f = (x(1)**2 - x(2))**2 + (1.0_dp + x(1))**2
ELSE IF (nprob == 5) THEN

//     Test problem 5 (Intermediate Rosenbrock)

  f = 10.0_dp*(x(1)**2 - x(n))**2 + (1.0_dp + x(1))**2
ELSE IF (nprob == 6) THEN

//     Test problem 6 (Equation (9.1.15) in Fletcher's book)

  f = - x(1) - x(2)
  con(1) = x(2) - x(1)**2
  con(2) = 1.0_dp - x(1)**2 - x(2)**2
ELSE IF (nprob == 7) THEN

//     Test problem 7 (Equation (14.4.2) in Fletcher's book)

  f = x(3)
  con(1) = 5.0_dp*x(1) - x(2) + x(3)
  con(2) = x(3) - x(1)**2 - x(2)**2 - 4.0_dp*x(2)
  con(m) = x(3) - 5.0_dp*x(1) - x(2)
ELSE IF (nprob == 8) THEN

//     Test problem 8 (Rosen-Suzuki)

  f = x(1)**2 + x(2)**2 + 2.0*x(3)**2 + x(4)**2 - 5.0_dp*x(1) - 5.0_dp*x(2)  &
      - 21.0_dp*x(3) + 7.0_dp*x(4)
  con(1) = 8.0_dp - x(1)**2 - x(2)**2 - x(3)**2 - x(4)**2 - x(1) + x(2)  &
           - x(3) + x(4)
  con(2) = 10._dp - x(1)**2 - 2._dp*x(2)**2 - x(3)**2 - 2._dp*x(4)**2 + x(1) + x(4)
  con(m) = 5.0_dp - 2.0*x(1)**2 - x(2)**2 - x(3)**2 - 2.0_dp*x(1) + x(2) + x(4)
ELSE IF (nprob == 9) THEN

//     Test problem 9 (Hock and Schittkowski 100)

  f = (x(1) - 10._dp)**2 + 5._dp*(x(2) - 12._dp)**2 + x(3)**4 + 3._dp*(x(4)  &
       - 11._dp)**2 + 10._dp*x(5)**6 + 7._dp*x(6)**2 + x(7)**4 - 4._dp*x(6)*x(7) &
       - 10._dp*x(6) - 8._dp*x(7)
  con(1) = 127._dp - 2._dp*x(1)**2 - 3._dp*x(2)**4 - x(3) - 4._dp*x(4)**2   &
           - 5._dp*x(5)
  con(2) = 282._dp - 7._dp*x(1) - 3._dp*x(2) - 10._dp*x(3)**2 - x(4) + x(5)
  con(3) = 196._dp - 23._dp*x(1) - x(2)**2 - 6._dp*x(6)**2 + 8._dp*x(7)
  con(4) = - 4._dp*x(1)**2 - x(2)**2 + 3._dp*x(1)*x(2) - 2._dp*x(3)**2 - 5._dp*x(6)  &
           + 11._dp*x(7)
ELSE IF (nprob == 10) THEN

//     Test problem 10 (Hexagon area)

  f = - 0.5_dp*(x(1)*x(4) - x(2)*x(3) + x(3)*x(n) - x(5)*x(n) + x(5)*x(8) &
                - x(6)*x(7))
  con(1) = 1.0_dp - x(3)**2 - x(4)**2
  con(2) = 1.0_dp - x(n)**2
  con(3) = 1.0_dp - x(5)**2 - x(6)**2
  con(4) = 1.0_dp - x(1)**2 - (x(2) - x(n))**2
  con(5) = 1.0_dp - (x(1) - x(5))**2 - (x(2) - x(6))**2
  con(6) = 1.0_dp - (x(1) - x(7))**2 - (x(2) - x(8))**2
  con(7) = 1.0_dp - (x(3) - x(5))**2 - (x(4) - x(6))**2
  con(8) = 1.0_dp - (x(3) - x(7))**2 - (x(4) - x(8))**2
  con(9) = 1.0_dp - x(7)**2 - (x(8) - x(n))**2
  con(10) = x(1)*x(4) - x(2)*x(3)
  con(11) = x(3)*x(n)
  con(12) = - x(5)*x(n)
  con(13) = x(5)*x(8) - x(6)*x(7)
  con(m) = x(n)
END IF

RETURN
END SUBROUTINE calcfc
*/
