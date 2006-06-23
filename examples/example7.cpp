#include <aa_aaf.h>

/* By using affine arithmetic we should be able to make an estimate of
 * the variance of a function over an input range. */

unsigned get_dominant_term(const AAF& x) {
    unsigned best = 0;
    double best_size = 0;
    for(unsigned i = 0; i < x.get_length(); i++) {
        double cf = fabs(x.get_coeff(i));
        if(cf > best_size) {
            best_size = cf;
            best = x.get_index(i);
        }
    }
    return best;
}

// Measure the derivative of result wrt x
double measure_dim_gradient(const AAF& result, const AAF& x) {
    unsigned dom = get_dominant_term(x);
    printf("dom: %d\n", dom);
    return result.index_coeff(dom) / x.index_coeff(dom);
}

// idea is to bound the interval around zero using the extents of result
AAF project(const AAF& result, const AAF& x) {
    return x - result / measure_dim_gradient(result, x);
}

int main() {
    AAF dummy1(interval(-1e-10, 1e-10)); // move away from zero to
                                         // help detect bugs.

    AAF x(interval(1.2, 1.5));

    AAF result = x*x - 2;

    std::cout << "x = " << x
              << std::endl;

    std::cout << result
              << (result.straddles_zero()?
                  " contains zero":
                  " does not contain zero")
              << std::endl;

    std::cout << "result_x = " << measure_dim_gradient(result, x)
              << std::endl;

    std::cout << "new x = " << project(result, x)
              << project(result, x).convert()
              << std::endl;

    std::cout << "Trial 2, 2 dimensions\n";

    AAF y(interval(0.5, 1.));

    result = x*y - 1;

    std::cout << "y = " << y
              << std::endl;

    std::cout << result
              << (result.straddles_zero()?
                  " contains zero":
                  " does not contain zero")
              << std::endl;

    std::cout << "result_x = " << measure_dim_gradient(result, x)
              << std::endl;
    std::cout << "result_y = " << measure_dim_gradient(result, y)
              << std::endl;

    AAF proj = project(project(result, x), y) * x;

    std::cout << "new x = " << proj << proj.convert()
              << std::endl;

    return 0;
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/


// vim: filetype=c++:expandtab:shiftwidth=4:tabstop=8:softtabstop=4 :
