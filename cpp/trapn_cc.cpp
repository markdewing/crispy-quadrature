
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>

const double pi = M_PI;

double transform_cc(double t)
{
    double L = 1.0;

    // Can cause performance problems (underflow?) near pi/2
    //return L/std::tan(t);

    // Use alternate formulation for cotangent
    double sin_t;
    double cos_t;
    sincos(t, &sin_t, &cos_t);
    return L*cos_t/sin_t;
}

double jacobian_cc(double t)
{
    double L = 1.0;
    double s = std::sin(t);
    return L/(s*s);
}

const int ndim = 6;

double fn(double xx[ndim])
{
    double x = xx[0];
    double y = xx[1];
    double z = xx[2];
    double x2 = xx[3];
    double y2 = xx[4];
    double z2 = xx[5];
    return exp(-x*x) * exp(-2*y*y) * exp(-z*z) * exp(-x2*x2)* exp(-y2*y2) * exp(-z2*z2);
}

// https://stackoverflow.com/questions/46782444/how-to-convert-a-linear-index-to-subscripts-with-support-for-negative-strides
void ind2sub(int idx, const int* shape, int* indices)
{
    for (int i = 0; i < ndim; i++)
    {
        int s = idx % shape[i];
        idx -= s;
        idx /= shape[i];
        indices[i] = s;
    }
}


template <typename F>
double trapn_cc(int n, F &f)
{
    double h = pi/(n+1);

    double total = 0.0;

    int npts = std::pow(n, ndim);
    int indices[ndim];
    int nn1[ndim];
    for (int i = 0; i < ndim; i++) {
        nn1[i] = n;
    }

#pragma omp parallel for reduction(+:total)
    for (int i = 0; i < npts; i++) {
        double xx[ndim];
        ind2sub(i, nn1, indices);
        double jac = 1.0;
        for (int j = 0; j < ndim; j++) {
            // Map to (0,pi)
            double x = (indices[j]+1) * h;

            // Map to (-oo, oo)
            xx[j] = transform_cc(x);

            jac *= jacobian_cc(x);
        }
        total += jac*f(xx);
    }

    return total*std::pow(h, ndim);
}

int main()
{
    //int n = 16;
    double expected = std::sqrt(pi/2) * std::pow(std::sqrt(pi), 5);
    std::cout << "# expected = " << expected << std::endl;
    std::cout << "# n      val     diff    time" << std::endl;
    for (int n = 10; n <= 24; n++) {
        double start = omp_get_wtime();
        double val = trapn_cc(n, fn);
        double end = omp_get_wtime();
        double diff = expected-val;
        std::cout << n << " " << val << " " << diff << " " << end-start << std::endl;
    }

}
