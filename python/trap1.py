
from math import exp

# Basic trapezoidal rule from a to b with n divisions

def trap1(a, b, n, f1):
    # Size of each interval
    h = (1.0/n)*(b-a)

    # first interval
    total = 1/2.0*h*f1(a)
    # last interval
    total += 1/2.0*h*f1(b)

    for i in range(1,n):
        total += h*f1(a+h*i)

    return total

# Integral is 1/2 x**2
def linear_f(x):
    return x

# Integral is 1/3 x**3
def quadratic_f(x):
    return x*x

# Integral is -exp(-x)
def exp_f(x):
    return exp(-x)


if __name__ == '__main__':
    a = 0.0
    b = 1.0
    n = 100

    # ----------------
    # Linear integrand
    # ----------------
    #expected_val = 0.5*(b*b -a*a)
    #val = trap1(a, b, n, linear_f)
    #print(val,(expected_val-val))

    # -------------------
    # Quadratic integrand
    # -------------------
    expected_val = (b**3 - a**3)/3
    val = trap1(a, b, n, quadratic_f)
    print("expected value =",val," diff from expected =",(expected_val-val))

    # ---------------------
    # Exponential integrand
    # ---------------------
    #expected_val = -(exp_f(b) - exp_f(a))
    #val = trap1(a, b, n, exp_f)
    #print(val,(expected_val-val))

