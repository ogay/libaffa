#include "aa.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

static int failures = 0;

static void require(bool ok, const char *message) {
    if (!ok) {
        std::cerr << message << '\n';
        ++failures;
    }
}

// Check each value against the affine form with the input noise fixed.
// This also checks preserved correlations, beyond the converted interval.
static void check(double a, double b) {
    const AAF input = AAF(interval(a, b));
    const AAF output = exp(input);
    const interval bounds = output.convert();
    const unsigned source = input.get_index(0);
    double error = 0;
    for (unsigned j = 0; j < output.get_length(); ++j)
        if (output.get_index(j) != source)
            error += std::fabs(output.get_coeff(j));
    for (int k = 0; k <= 2048; ++k) {
        const long double e = -1.L + 2.L*k/2048;
        const long double x = input.get_center() + e*input.get_coeff(0);
        const long double value = std::exp(x);
        const long double line = output.get_center() + e*output.index_coeff(source);
        // The library does not track rounding errors: this test targets the
        // approximation defect, allowing a small floating-point tolerance.
        const long double tol = 2e-12L*std::max(1.L, std::fabs(value));
        if (!(bounds.left() <= value + tol && bounds.right() >= value - tol
              && std::fabs(value-line) <= error + tol)) {
            std::cerr << "Failure on [" << a << ", " << b << "] at " << x
                      << ": range [" << bounds.left() << ", " << bounds.right()
                      << "], value " << value << ", affine residual "
                      << std::fabs(value-line) << ", error " << error << '\n';
            ++failures;
            return;
        }
    }
}

int main() {
    check(0, 1);
    check(-4, -1);
    check(1, 4);
    check(-2, 3);
    check(0, 0);
    check(1, 1);
    check(-800, -799); // both endpoints underflow
    check(0, 1e-17);  // endpoint exponentials round to the same value
    for (int i = -20; i <= 20; ++i)
        check(i/4.0, i/4.0 + 0.125);
    require(exp(AAF(0.)).get_center() == 1, "exp(0) must equal one");
    return failures ? 1 : 0;
}
