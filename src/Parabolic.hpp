#pragma once

#include <cmath>

#include "Constants.hpp"
#include "Solver.hpp"
#include "Vector.hpp"

class Parabolic {
public:
    /// @param wavenumber   волновое число k = 2*pi/lambda
    /// @param amplitude    амплитуда падающей волны E0
    /// @param focalLength  фокусное расстояние f0
    /// @param mirrorRadius радиус зеркала r_m (определяет верхний предел по tau)
    /// @param initialPhase начальная фаза psi0
    Parabolic(double wavenumber, double amplitude, double focalLength, double mirrorRadius, double initialPhase)
        : k_(wavenumber),
          w_(wavenumber * constants::kSpeedOfLight),
          psi0_(initialPhase),
          e0_(amplitude),
          f0_(focalLength),
          rm_(mirrorRadius),
          sigma_(std::sqrt(2.0 * focalLength)) {}

    void setConfig(const Vector3& obsPoint, double time) {
        xp_ = obsPoint.x;
        yp_ = obsPoint.y;
        zp_ = obsPoint.z;
        t_ = time;
    }

    Vector3 integrateE(const Vector3& obsPoint, double time, int n, int m) {
        setConfig(obsPoint, time);
        return Solver().integrateVectorField(
            0.0, tauMax(), zeroLimit(), fullCircleLimit(),
            [this](double tau, double fi) { return computeFe(tau, fi); }, n, m);
    }

    Vector3 integrateH(const Vector3& obsPoint, double time, int n, int m) {
        setConfig(obsPoint, time);
        return Solver().integrateVectorField(
            0.0, tauMax(), zeroLimit(), fullCircleLimit(),
            [this](double tau, double fi) { return computeFh(tau, fi); }, n, m);
    }

    Vector3 computeBe(double tau, double fi) const {
        const double x = mirrorX(tau);
        const double y = mirrorY(tau, fi);
        const double z = mirrorZ(tau, fi);
        return Vector3(-2.0 * y * (x - xp_),
                       -2.0 * y * (y - yp_),
                       -2.0 * y * (z - zp_));
    }

    Vector3 computeA(double tau, double fi) const {
        return Vector3(2.0 * mirrorY(tau, fi), 4.0 * f0_, 0.0);
    }

    Vector3 computeBh(double tau, double fi) const {
        const double x = mirrorX(tau);
        const double y = mirrorY(tau, fi);
        const double z = mirrorZ(tau, fi);
        return Vector3(4.0 * f0_ * (z - zp_),
                       -2.0 * y * (z - zp_),
                       2.0 * y * (y - yp_) - 4.0 * f0_ * (x - xp_));
    }

    double computeU(double tau, double fi) const {
        const double dx = mirrorX(tau) - xp_;
        const double dy = mirrorY(tau, fi) - yp_;
        const double dz = mirrorZ(tau, fi) - zp_;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    Vector3 computeFe(double tau, double fi) const {
        const Vector3 be = computeBe(tau, fi);
        const Vector3 a = computeA(tau, fi);
        const double u = computeU(tau, fi);
        const double arg = k_ * (mirrorX(tau) - u) + w_ * t_ + psi0_;
        const Vector3 sinpart = a * (k_ / u) + be * (k_ / (u * u));
        const Vector3 cospart = -be / (u * u * u);

        return (sinpart * std::sin(arg) + cospart * std::cos(arg)) *
               tau * e0_ / (4 * constants::kPi);
    }

    Vector3 computeFh(double tau, double fi) const {
        const Vector3 bh = computeBh(tau, fi);
        const double u = computeU(tau, fi);
        const double arg = k_ * (mirrorX(tau) - u) + w_ * t_ + psi0_;
        const Vector3 sinpart = bh * k_ / (u * u);
        const Vector3 cospart = -bh / (u * u * u);

        return (sinpart * std::sin(arg) + cospart * std::cos(arg)) *
               tau * e0_ / (4 * constants::kPi);
    }

private:
    double tauMax() const { return rm_ / std::sqrt(2 * f0_); }

    static Solver::Limit zeroLimit() {
        return [](double) { return 0.0; };
    }
    static Solver::Limit fullCircleLimit() {
        return [](double) { return 2.0 * constants::kPi; };
    }

    double mirrorX(double tau) const {
        return (tau * tau - sigma_ * sigma_) / 2;
    }
    double mirrorY(double tau, double fi) const {
        return sigma_ * tau * std::cos(fi);
    }
    double mirrorZ(double tau, double fi) const {
        return sigma_ * tau * std::sin(fi);
    }

    double k_;      ///< волновое число
    double w_;      ///< циклическая частота
    double psi0_;   ///< начальная фаза
    double e0_;     ///< амплитуда
    double f0_;     ///< фокусное расстояние
    double rm_;     ///< радиус зеркала
    double sigma_;  ///< sqrt(2 f0)
    double t_ = 0.0;   ///< время

    double xp_ = 0.0;  ///< точка наблюдения
    double yp_ = 0.0;
    double zp_ = 0.0;
};