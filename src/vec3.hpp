#pragma once
#include <array>
#include <cmath>

template<typename T>
using Vec3 = std::array<T, 3>;

Vec3<double> operator*(double scalar, Vec3<double> x)
{
	for (auto& v : x)
		v *= scalar;
	return x;
}

Vec3<double> operator+(Vec3<double> x, const Vec3<double>& y)
{
	for (std::size_t i = 0; i < x.size(); ++i)
		x[i] += y[i];
	return x;
}

double operator*(const Vec3<double>& x, const Vec3<double>& y)
{
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

Vec3<double> operator^(const Vec3<double>& x, const Vec3<double>& y)
{
	return {x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]};
}

double norm_sq(const Vec3<double>& x)
{
	return x * x;
}

double norm(const Vec3<double>& x)
{
	return std::sqrt(norm_sq(x));
}
