#pragma once

/// Трёхмерный вектор с базовой векторной алгеброй.
struct Vector3 {
    double x, y, z;

    explicit Vector3(double x = 0.0, double y = 0.0, double z = 0.0)
        : x(x), y(y), z(z) {}

    Vector3 operator*(double a) const { return Vector3(x * a, y * a, z * a); }

    Vector3 operator/(double a) const { return (*this) * (1 / a); }

    Vector3 operator+(const Vector3& a) const {
        return Vector3(x + a.x, y + a.y, z + a.z);
    }

    Vector3 operator-() const { return (*this) * (-1); }

    Vector3 operator-(const Vector3& a) const { return (*this) + (-a); }
};
