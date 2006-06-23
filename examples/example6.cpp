#include "aa.h"

/* A toy example from Easyval library to demonstrate underflow.
 * Ideally we would be able to detect this error. */

const double eps = 1e-10;

int main() {
  AAF a(interval(77617-eps, 77617+eps));
  AAF b(interval(33096-eps, 33096+eps));

  AAF result =  333.75*pow( b, 6 )
               + a*a*(11*a*a*b*b - pow( b, 6 ) - 121*pow( b, 4 ) - 2)
               + 5.5*pow( b, 8 ) + a/(2*b);
  std::cout << result << std::endl;

  std::cout << "Although the input values("
        << a.get_center() << ", " << b.get_center()
        << ") only had an uncertainty of " << eps
        << ", the output interval is:" << result.convert()
        << ".  The correct value is actually about -0.82740."
        << std::endl;
  return 0;
}
