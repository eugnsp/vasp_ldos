#pragma once
#include <array>
#include <cmath>
#include <ostream>

template<typename T>
using Vec3 = std::array<T, 3>;

template<typename T>
using Basis3 = std::array<Vec3<T>, 3>;

Vec3<double> operator*(double scalar, Vec3<double> vec)
{
	for (auto& v : vec)
		v *= scalar;
	return vec;
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

double norm_sq(const Vec3<double>& vec)
{
	return vec * vec;
}

double norm(const Vec3<double>& vec)
{
	return std::sqrt(norm_sq(vec));
}

template<class T>
std::ostream& operator<<(std::ostream& os, const Vec3<T>& vec)
{
	os << '(' << vec[0] << ", " << vec[1] << ", " << vec[2] << ')';
	return os;
}
