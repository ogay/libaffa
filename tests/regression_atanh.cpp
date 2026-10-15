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
    const AAF output = atanh(input);
    const interval bounds = output.convert();
    const unsigned source = input.get_index(0);
    double error = 0;
    for (unsigned j = 0; j < output.get_length(); ++j)
        if (output.get_index(j) != source)
            error += std::fabs(output.get_coeff(j));
    for (int k = 0; k <= 2048; ++k) {
        const long double e = -1.L + 2.L*k/2048;
        const long double x = input.get_center() + e*input.get_coeff(0);
        const long double value = std::atanh(x);
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
    check(-0.5, 0.5);
    check(0.1, 0.9);
    check(-0.9, -0.1);
    check(-0.99, 0.99);
    check(0.99, 0.999);
    check(0, 0);
    check(0.5, 0.5);
    check(-0.5, -0.5);
    for (int i = -90; i <= 80; i += 5)
        check(i/100.0, (i+10)/100.0);
    require(atanh(AAF(0.)).get_center() == 0, "atanh(0) must equal zero");
    require(atanh(AAF(interval(-1, 0))).is_nan(), "reject left pole");
    require(atanh(AAF(interval(0, 1))).is_nan(), "reject right pole");
    require(atanh(AAF(interval(2, 3))).is_nan(), "reject invalid domain");
    require(atanh(AAF(AAF_TYPE_NAN)).is_nan(), "NaN must propagate");
    require(atanh(AAF(AAF_TYPE_INFINITE)).is_nan(), "reject unbounded domain");
    return failures ? 1 : 0;
}
