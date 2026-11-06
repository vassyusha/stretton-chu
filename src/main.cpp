#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "Cartesian.hpp"
#include "Constants.hpp"
#include "Parabolic.hpp"

namespace {

void printVector(const char* name, const Vector3& v) {
    std::cout << "  " << name << " = (" << v.x << ", " << v.y << ", " << v.z
              << ")\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    using constants::kPi;
    using constants::kSpeedOfLight;

    const int n = (argc > 1) ? std::atoi(argv[1]) : 1000;

    // Параметры задачи (СГС, длины в см).
    const double rm = 0.01;                  // радиус зеркала
    const double f0 = 2.0 * rm;              // фокусное расстояние (f# = 1)
    const double lambda = 1e-4;              // длина волны
    const double k = 2 * kPi / lambda;       // волновое число
    const double psi0 = kPi / 2.0;           // начальная фаза
    const double power = 1e22;               // см. формулу E0 ниже
    const double e0 = std::sqrt(2.0 * power / (kSpeedOfLight * rm * rm));
    const double time = 0.0;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "E0 = " << e0 << ", grid = " << n << " x " << n << "\n";

    Cartesian cartesian(k, e0, f0, rm, psi0);
    Parabolic parabolic(k, e0, f0, rm, psi0);

    const Vector3 points[] = {
        Vector3(0.0, 0.0, 0.0),
        Vector3(0.0, 0.5 * lambda, 0.5 * lambda),
    };

    for (const Vector3& p : points) {
        std::cout << "\nobservation point (" << p.x << ", " << p.y << ", "
                  << p.z << ")\n";

        std::cout << "Cartesian:\n";
        printVector("E", cartesian.integrateE(p, time, n, n));
        printVector("H", cartesian.integrateH(p, time, n, n));

        std::cout << "Parabolic:\n";
        printVector("E", parabolic.integrateE(p, time, n, n));
        printVector("H", parabolic.integrateH(p, time, n, n));
    }

    return 0;
}