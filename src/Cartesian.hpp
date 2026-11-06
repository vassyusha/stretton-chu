#pragma once

#include <cmath>

#include "Constants.hpp"
#include "Solver.hpp"
#include "Vector.hpp"

class Cartesian {
public:
    /// @param wavenumber   волновое число k = 2*pi/lambda
    /// @param amplitude    амплитуда падающей волны E0
    /// @param focalLength  фокусное расстояние f0
    /// @param mirrorRadius радиус зеркала r_m (радиус области интегрирования)
    /// @param initialPhase начальная фаза psi0
    Cartesian(double wavenumber, double amplitude, double focalLength, double mirrorRadius, double initialPhase)
        : k_(wavenumber),
          w_(wavenumber * constants::kSpeedOfLight),
          psi0_(initialPhase),
          e0_(amplitude),
          f0_(focalLength),
          rm_(mirrorRadius) {}

    void setConfig(const Vector3& obsPoint, double time) {
        xp_ = obsPoint.x;
        yp_ = obsPoint.y;
        zp_ = obsPoint.z;
        t_ = time;
    }

    Vector3 integrateE(const Vector3& obsPoint, double time, int n, int m) {
        setConfig(obsPoint, time);
        return Solver().integrateVectorField(
            -rm_, rm_, lowerZ(), upperZ(),
            [this](double y, double z) { return computeFe(y, z); }, n, m);
    }

    Vector3 integrateH(const Vector3& obsPoint, double time, int n, int m) {
        setConfig(obsPoint, time);
        return Solver().integrateVectorField(
            -rm_, rm_, lowerZ(), upperZ(),
            [this](double y, double z) { return computeFh(y, z); }, n, m);
    }

    Vector3 computeBe(double y, double z) const {
        const double x = mirrorX(y, z);
        return Vector3(-2.0 * y * (x - xp_),
                       -2.0 * y * (y - yp_),
                       -2.0 * y * (z - zp_));
    }

    Vector3 computeA(double y, double /*z*/) const {
        return Vector3(2.0 * y, 4.0 * f0_, 0.0);
    }

    Vector3 computeBh(double y, double z) const {
        const double x = mirrorX(y, z);
        return Vector3(4.0 * f0_ * (z - zp_),
                       -2.0 * y * (z - zp_),
                       2.0 * y * (y - yp_) - 4.0 * f0_ * (x - xp_));
    }

    double computeU(double y, double z) const {
        const double dx = mirrorX(y, z) - xp_;
        const double dy = y - yp_;
        const double dz = z - zp_;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    Vector3 computeFe(double y, double z) const {
        const Vector3 a = computeA(y, z);
        const Vector3 be = computeBe(y, z);
        const double u = computeU(y, z);
        const double arg = k_ * (mirrorX(y, z) - u) + w_ * t_ + psi0_;
        const Vector3 sinpart = a * (k_ / u) + be * (k_ / (u * u));
        const Vector3 cospart = -be / (u * u * u);

        return (sinpart * std::sin(arg) + cospart * std::cos(arg)) /
               (2.0 * f0_) * e0_ / (4 * constants::kPi);
    }

    Vector3 computeFh(double y, double z) const {
        const Vector3 bh = computeBh(y, z);
        const double u = computeU(y, z);
        const double arg = k_ * (mirrorX(y, z) - u) + w_ * t_ + psi0_;
        const Vector3 sinpart = bh * k_ / (u * u);
        const Vector3 cospart = -bh / (u * u * u);

        return (sinpart * std::sin(arg) + cospart * std::cos(arg)) /
               (2.0 * f0_) * e0_ / (4 * constants::kPi);
    }

private:
    Solver::Limit lowerZ() const {
        const double rm = rm_;
        return [rm](double y) { return -std::sqrt(rm * rm - y * y); };
    }
    Solver::Limit upperZ() const {
        const double rm = rm_;
        return [rm](double y) { return std::sqrt(rm * rm - y * y); };
    }

    double mirrorX(double y, double z) const {
        return (y * y + z * z) / (4.0 * f0_) - f0_;
    }

    double k_;     ///< волновое число
    double w_;     ///< циклическая частота
    double psi0_;  ///< начальная фаза
    double e0_;    ///< амплитуда
    double f0_;    ///< фокусное расстояние
    double rm_;    ///< радиус зеркала
    double t_ = 0.0;   ///< время

    double xp_ = 0.0;  ///< точка наблюдения
    double yp_ = 0.0;
    double zp_ = 0.0;
};