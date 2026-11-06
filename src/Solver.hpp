#pragma once

#include <functional>

#include "Vector.hpp"

/// Численное интегрирование функций двух переменных по области вида
///   a <= x <= b,  y1(x) <= y <= y2(x).
class Solver {
public:
    using Limit = std::function<double(double)>;
    using Scalar2D = std::function<double(double, double)>;
    using Vector2D = std::function<Vector3(double, double)>;

    /// Метод трапеций (n узлов по первой переменной, m по второй).
    double trapezoidRule(double lower1, double upper1,
                         const Limit& lower2, const Limit& upper2,
                         const Scalar2D& func, int n, int m) const;

    /// Метод ячеек: ячейка учитывается, если её центр лежит в области.
    double cellMethodCent(double lower1, double upper1,
                          const Limit& lower2, const Limit& upper2,
                          const Scalar2D& func, int n, int m) const;

    /// Метод ячеек: ячейка учитывается, если все её углы лежат в области.
    double cellMethodCorn(double lower1, double upper1,
                          const Limit& lower2, const Limit& upper2,
                          const Scalar2D& func, int n, int m) const;

    /// Покомпонентное интегрирование векторного поля методом трапеций.
    Vector3 integrateVectorField(double lower1, double upper1,
                                 const Limit& lower2, const Limit& upper2,
                                 const Vector2D& vecFunc, int n, int m) const;

private:
    static bool isInCorners(double lower1, double upper1,
                            const Limit& lower2, const Limit& upper2,
                            double minLimit2, int i, int j,
                            double dx1, double dx2);

    static bool isInCentre(double lower1, double upper1,
                           const Limit& lower2, const Limit& upper2,
                           double x, double y);
};
