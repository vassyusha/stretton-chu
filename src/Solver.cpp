#include "Solver.hpp"

#include <algorithm>
#include <limits>

namespace {

/// Интеграл по второй переменной при фиксированной первой (формула трапеций).
double lineIntegral(double x, const Solver::Limit& lower2,
                    const Solver::Limit& upper2,
                    const Solver::Scalar2D& func, int m) {
    const double lo = lower2(x);
    const double hi = upper2(x);
    const double dx2 = (hi - lo) / (m - 1);

    double sum = (func(x, lo) + func(x, hi)) / 2;
    for (int j = 1; j < m - 1; j++) {
        sum += func(x, lo + j * dx2);
    }
    return sum * dx2;
}

}  // namespace

double Solver::trapezoidRule(double lower1, double upper1,
                             const Limit& lower2, const Limit& upper2,
                             const Scalar2D& func, int n, int m) const {
    const double dx1 = (upper1 - lower1) / (n - 1);

    // Граничные линии x = lower1 и x = upper1 входят с весом 1/2.
    double integral = lineIntegral(lower1, lower2, upper2, func, m);
    integral += lineIntegral(upper1, lower2, upper2, func, m);
    integral /= 2;

    for (int i = 1; i < n - 1; i++) {
        integral += lineIntegral(lower1 + i * dx1, lower2, upper2, func, m);
    }
    return integral * dx1;
}

double Solver::cellMethodCent(double lower1, double upper1,
                              const Limit& lower2, const Limit& upper2,
                              const Scalar2D& func, int n, int m) const {
    const double dx1 = (upper1 - lower1) / (n - 1);

    double maxLimit2 = -std::numeric_limits<double>::infinity();
    double minLimit2 = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
        maxLimit2 = std::max(maxLimit2, upper2(lower1 + i * dx1));
        minLimit2 = std::min(minLimit2, lower2(lower1 + i * dx1));
    }
    const double dx2 = (maxLimit2 - minLimit2) / (m - 1);

    double integral = 0.0;
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < m - 1; j++) {
            const double xav = lower1 + (i + 1.0 / 2) * dx1;
            const double yav = minLimit2 + (j + 1.0 / 2) * dx2;
            if (isInCentre(lower1, upper1, lower2, upper2, xav, yav)) {
                integral += func(xav, yav);
            }
        }
    }
    return integral * dx1 * dx2;
}

double Solver::cellMethodCorn(double lower1, double upper1,
                              const Limit& lower2, const Limit& upper2,
                              const Scalar2D& func, int n, int m) const {
    const double dx1 = (upper1 - lower1) / (n - 1);

    double maxLimit2 = -std::numeric_limits<double>::infinity();
    double minLimit2 = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
        maxLimit2 = std::max(maxLimit2, upper2(lower1 + i * dx1));
        minLimit2 = std::min(minLimit2, lower2(lower1 + i * dx1));
    }
    const double dx2 = (maxLimit2 - minLimit2) / (m - 1);

    double integral = 0.0;
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < m - 1; j++) {
            if (isInCorners(lower1, upper1, lower2, upper2, minLimit2,
                            i, j, dx1, dx2)) {
                const double xav = lower1 + (i + 1.0 / 2) * dx1;
                const double yav = minLimit2 + (j + 1.0 / 2) * dx2;
                integral += func(xav, yav);
            }
        }
    }
    return integral * dx1 * dx2;
}

bool Solver::isInCorners(double lower1, double upper1,
                         const Limit& lower2, const Limit& upper2,
                         double minLimit2, int i, int j,
                         double dx1, double dx2) {
    const double x0 = lower1 + i * dx1;
    const double x1 = lower1 + (i + 1) * dx1;

    const double y0 = minLimit2 + j * dx2;
    const double y1 = minLimit2 + (j + 1) * dx2;

    return x0 >= lower1 && x1 <= upper1 &&
           y0 >= lower2(x0) && y0 >= lower2(x1) &&
           y1 <= upper2(x0) && y1 <= upper2(x1);
}

bool Solver::isInCentre(double lower1, double upper1,
                        const Limit& lower2, const Limit& upper2,
                        double x, double y) {
    return x > lower1 && x < upper1 && y > lower2(x) && y < upper2(x);
}

Vector3 Solver::integrateVectorField(double lower1, double upper1,
                                     const Limit& lower2, const Limit& upper2,
                                     const Vector2D& vecFunc,
                                     int n, int m) const {
    Vector3 result(0.0, 0.0, 0.0);

    // Каждая компонента интегрируется отдельно.
    result.x = trapezoidRule(lower1, upper1, lower2, upper2,
        [&vecFunc](double a, double b) { return vecFunc(a, b).x; }, n, m);
    result.y = trapezoidRule(lower1, upper1, lower2, upper2,
        [&vecFunc](double a, double b) { return vecFunc(a, b).y; }, n, m);
    result.z = trapezoidRule(lower1, upper1, lower2, upper2,
        [&vecFunc](double a, double b) { return vecFunc(a, b).z; }, n, m);

    return result;
}
