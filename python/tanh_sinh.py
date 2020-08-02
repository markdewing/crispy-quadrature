
from math import exp,pi,tanh,sinh,cosh,log,sqrt
from scipy.special import lambertw

# From https://arxiv.org/abs/2007.15057
# Tanh-sinh quadrature for single and multiple integration
#   using floating-point arithmetic
#  Joren Vanherck, Bart Soree, Wim Magnus


# Trapezoidal rule from a to b with n divisions
# For comparison

def trap_base(a, b, n, f1):

    # Size of each interval
    h = (1.0/n)*(b-a)

    # first interval
    total = 1/2.0*h*f1(a)
    # last interval
    total += 1/2.0*h*f1(b)

    for i in range(1,n):
        x = a + h*i
        total += h*f1(x)

    return total


# Trapezoidal rule for use with transformed coordinates

def trap1(h, n, f1, tr, wt):

    total = h*wt(0.0)*f1(tr(0.0))
    for i in range(1,n):
        x = h*i
        total += h*f1(tr(x))*wt(x)
        x = -h*i
        total += h*f1(tr(x))*wt(x)

        # Spacing diagnostics
        #yi_p = exp(0.5*pi*sinh(x))/cosh(0.5*pi*sinh(x))
        #yi_m = exp(-0.5*pi*sinh(x))/cosh(0.5*pi*sinh(x))
        #print('x = ',x,yi_p,yi_m)

    return total

def psi(t):
    return tanh(0.5*pi*sinh(t))

def psi_prime(t):
    return 0.5*pi*cosh(t)/cosh(0.5*pi*sinh(t))**2

def tanh_quad(n,f1):
   # Optimal h regardless of FP precision
   #h = 2*lambertw(pi*n).real/n

   # Maximum value to use to ensure spacing
   #  at the endpoints is large enough
   # From the paper, for double precision.
   tmax = 6.112
   h = tmax/n
   print('optimal h = ',h)
   return trap1(h,n,f1,psi,psi_prime)
    

# Some test integrands

# Integral is 1/2 x**2
def linear_f(x):
    return x

# Integral is 1/3 x**3
def quadratic_f(x):
    return x*x

# Integral is -exp(-x)
def exp_f(x):
    return exp(-x)

def exp_sq(x):
    return exp(-x*x)

def inv_sqrt(x):
    return 1.0/sqrt(x)

def inv_abs(x):
    return 1.0/abs(x)

def transform(a,b,f):
    def ft(x):
        y = 0.5*((b-a)*x + (b+a))
        return f(y)
    return ft


if __name__ == '__main__':

    a =  0.0 
    b =  10.0
    jac = (b-a)*0.5

    n = 20

    # ----------------
    # Linear integrand
    # ----------------
    #expected_val = 0.5*(b*b -a*a)
    #val = trap_base(a, b, n, linear_f)
    #print(val,(expected_val-val))
    #linear_tr = transform(a,b,linear_f)
    #val = jac*tanh_quad(n, linear_tr)
    #print(val,(expected_val-val))

    # -------------------
    # Quadratic integrand
    # -------------------
    #expected_val = (b**3 - a**3)/3
    #val = trap_base(a, b, n, quadratic_f)
    #print(val,(expected_val-val))
    #val = tanh_quad(n, quadratic_f)
    #print(val,(expected_val-val))



    # ----------------------
    # Exponential integrand
    # ----------------------

    expected_val = -(exp_f(b) - exp_f(a))

    val = trap_base(a, b, n, exp_f)
    print('trapezoidal rule value = ',val,' diff from expected value = ',(expected_val-val))

    exp_tr = transform(a,b,exp_f)
    val = jac*tanh_quad(n, exp_tr)
    print('tanh-sinh rule value = ',val, ' diff from expected value = ',(expected_val-val))

    # -----------------------------
    # Exponential squared integrand
    # -----------------------------

    # For infinite integration limits
    #expected_val = sqrt(pi)

    #val1 = trap_base(a, b, n, exp_sq)
    #print(val1,(expected_val-val1))

    #exp_sq_tr = transform(a,b,exp_sq)
    #val2 = jac*tanh_quad(n, exp_sq_tr)
    #print(val2,(expected_val-val2))

    # -------------------
    # Inverse square root
    # -------------------
    #inv_tr = transform(a,b,inv_abs)

    #val1 = jac*trap_base(a, b, n, inv_tr)
    #val2 = jac*tanh_quad(n, inv_tr)
    #print(val1,val2,(val1-val2))

