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
    const AAF output = sin(input);
    const interval bounds = output.convert();
    const unsigned source = input.get_index(0);
    double error = 0;
    for (unsigned j = 0; j < output.get_length(); ++j)
        if (output.get_index(j) != source)
            error += std::fabs(output.get_coeff(j));
    for (int k = 0; k <= 2048; ++k) {
        const long double e = -1.L + 2.L*k/2048;
        const long double x = input.get_center() + e*input.get_coeff(0);
        const long double value = std::sin(x);
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
    check(-4.41, 1.21);
    check(-1, 1);
    check(0, 6.2);
    check(-10, 10);
    check(0, 0);
    check(2, 2);
    check(0, 1e-8);
    for (int i = -40; i <= 40; ++i) {
        check(i/4.0, i/4.0 + 0.125);
        check(i/4.0, i/4.0 + 5.62);
    }
    const interval s = sin(AAF(interval(-4.41, 1.21))).convert();
    require(s.left() <= -1, "the internal sine minimum must be enclosed");
    require(sin(AAF(AAF_TYPE_NAN)).is_nan(), "NaN must propagate");
    return failures ? 1 : 0;
}
