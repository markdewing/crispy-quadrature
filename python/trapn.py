
from math import exp, pi, tan, sin, sqrt
import numpy as np
# Trapezoidal rule from -oo to oo with n divisions
# https://en.wikipedia.org/wiki/Numerical_integration#Integrals_over_infinite_intervals

# In ndim dimensions

# Coordinate transformation
@np.vectorize
def transform_inf(t):
    return t/(1.0 - t*t)

# Jacobian of the transformation
@np.vectorize
def jacobian_inf(t):
    return (1.0 + t*t) / (1.0 - t*t)**2

# https://stackoverflow.com/questions/46782444/how-to-convert-a-linear-index-to-subscripts-with-support-for-negative-strides
def ind2sub(idx, shape, indices):
    for i in range(len(shape)):
        s = idx % shape[i]
        idx -= s
        idx //= shape[i]
        indices[i] = s
    return indices

# Assume the integrand goes to zero faster than 1/x^2, so the endpoints evaluate to zero.

# n is the number of interior grid points.
# The number of intervals is n+1
# Including the end points, the number of grid points is n+2
#   Since they contribute zero to the integral, we can skip them and
#   only use the interior grid points.

def trapn(ndim, n, fn):
    a = -1.0
    b = 1.0

    # Size of each interval
    h = 2.0 / (n+1)
    hval = h**ndim

    total = 0.0
    nn1 = [n]*ndim
    idx = np.zeros(ndim, dtype=np.int64)
    x = np.zeros(ndim)
    xt = np.zeros(ndim)

    npts = n**ndim
    for i in range(npts):
        idx = ind2sub(i, nn1, idx) + 1

        # Map to (-1.0, 1.0)
        x[:] = idx*h - 1.0

        # Map to (-oo,oo)
        xt[:] = transform_inf(x)

        jac = np.prod(jacobian_inf(x))
        total += jac*fn(xt)

    return total*hval

# Integral is -exp(-x)
def exp_fn(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    x5 = x[4]
    x6 = x[5]
    return exp(-x1*x1)*exp(-2*x2*x2)*exp(-3*x3*x3)*exp(-x4*x4)*exp(-x5*x5)*exp(-x6*x6)


if __name__ == '__main__':

    # Number of grid points used in each dimension
    n = 8

    e = sqrt(pi) * sqrt(pi/2) * sqrt(pi/3) * sqrt(pi) * sqrt(pi) * sqrt(pi)
    print('Expected value: ',e)
    # ---------------------
    # Exponential integrand
    # ---------------------
    ndim = 6
    val = trapn(ndim, n, exp_fn)
    print('Integrated value: ',val, ' Difference: ',val-e)

