
from math import exp, pi, tan, sin, sqrt

# Clenshaw-Curtis - trapezoidal rule from -oo to oo with n divisions
# https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature

# Coordinate transformation
def transform_cc(t):
    L = 1.0
    return L/tan(t)

# Jacobian of the transformation
def jacobian_cc(t):
    L = 1.0
    return L/sin(t)**2

def trap_cc(n, f1):
    a = 0.0
    b = pi
    # Size of each interval
    h = (1.0/n)*(pi)

    total = 0.0

    for i in range(1,n):
        x =  h*i
        total += h*f1(transform_cc(x))*jacobian_cc(x)

    return total

# Integral is -exp(-x)
def exp_f(x):
    return exp(-x*x)


if __name__ == '__main__':
    a = 0.0
    b = pi
    n = 40


    print('Expected value: ',sqrt(pi))
    # ---------------------
    # Exponential integrand
    # ---------------------
    val = trap_cc(n, exp_f)
    print('Integrated value: ',val, ' Difference: ',val-sqrt(pi))

