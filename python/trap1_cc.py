
from math import exp, pi, tan, sin, sqrt

# Clenshaw-Curtis - trapezoidal rule from -oo to oo with n points
# https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature

# Coordinate transformation
def transform_cc(t):
    # Problem-dependent convergence parameter
    L = 1.0
    return L/tan(t)

# Jacobian of the transformation
def jacobian_cc(t):
    # Problem-dependent convergence parameter
    L = 1.0
    return L/sin(t)**2

def trap_cc(n, f1):
    a = 0.0
    b = pi
    # Size of each interval
    h = (1.0/(n+1))*(pi)

    total = 0.0

    for i in range(1, n+1):
        # grid points in (0,pi)
        x =  h*i
        total += h*f1(transform_cc(x))*jacobian_cc(x)

    return total

# Integral is -exp(-x)
def exp_f(x):
    return exp(-x*x)


if __name__ == '__main__':
    a = 0.0
    b = pi


    # We assume the end points (integrand at infinity) evaluate to zero,
    # so only the interior points are used.

    # Number of interior grid points.
    n = 10

    print('Expected value: ',sqrt(pi))
    # ---------------------
    # Exponential integrand
    # ---------------------
    val = trap_cc(n, exp_f)
    print('Integrated value: ',val, ' Difference: ',val-sqrt(pi))

