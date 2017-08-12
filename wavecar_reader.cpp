#include "wavecar_reader.hpp"
#include <cmath>
#include <type_traits>
#include <stdexcept>
#include <cassert>

namespace
{
constexpr double PI = 3.141592653589793238463;

// The value of 2 * m / hbar^2 in units of 1 / (eV * Ang^2)
constexpr double TWO_M_OVER_HBAR_SQ = 0.262465831;

// Simple linear algebra routines
double scalar(const Vec3<double>& x, const Vec3<double>& y)
{
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

double norm_sq(const Vec3<double>& x)
{
	return scalar(x, x);
}

double norm(const Vec3<double>& x)
{
	return std::sqrt(norm_sq(x));
}

Vec3<double> cross(const Vec3<double>& x, const Vec3<double>& y)
{
	Vec3<double> r;
	r[0] = x[1] * y[2] - x[2] * y[1];
	r[1] = x[2] * y[0] - x[0] * y[2];
	r[2] = x[0] * y[1] - x[1] * y[0];

	return r;
}

double det(const Vec3<double>& a1, const Vec3<double>& a2, const Vec3<double>& a3)
{
	return scalar(a1, cross(a2, a3));
}

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

template<typename T>
typename std::make_signed<T>::type index_shift(T i, T i_max)
{
	using T_s = typename std::make_signed<T>::type;

	const auto i_max_s = static_cast<T_s>(i_max);
	auto index = static_cast<T_s>(i);
	if (index > i_max_s)
		index -= (2 * i_max_s + 1);

	return index;
}
}

Wavecar_reader::Wavecar_reader(const std::string& filename)
	: File_reader(filename, File_type::BINARY)
{
	file_.exceptions(std::ifstream::badbit | std::ifstream::failbit);

	read_header();
	compute_reciprocal();
}

bool Wavecar_reader::is_single_precision() const
{
	return precision_ == Precision::SINGLE;
}

bool Wavecar_reader::is_double_precision() const
{
	return precision_ == Precision::DOUBLE;
}

std::size_t Wavecar_reader::n_spins() const
{
	return n_spins_;
}

std::size_t Wavecar_reader::n_kpoints() const
{
	return n_kpoints_;
}

std::size_t Wavecar_reader::n_bands() const
{
	return n_bands_;
}

double Wavecar_reader::e_cut() const
{
	return e_cut_;
}

const Vec3<double>& Wavecar_reader::a0() const
{
	return a0_;
}

const Vec3<double>& Wavecar_reader::a1() const
{
	return a1_;
}

const Vec3<double>& Wavecar_reader::a2() const
{
	return a2_;
}

double Wavecar_reader::a0_norm() const
{
	return norm(a0_);
}

double Wavecar_reader::a1_norm() const
{
	return norm(a1_);
}

double Wavecar_reader::a2_norm() const
{
	return norm(a2_);
}

std::size_t Wavecar_reader::size_g0() const
{
	return 2 * max_g0_ + 1;
}

std::size_t Wavecar_reader::size_g1() const
{
	return 2 * max_g1_ + 1;
}

std::size_t Wavecar_reader::size_g2() const
{
	return 2 * max_g2_ + 1;
}

template<typename T>
void Wavecar_reader::get_kpoint_data(
	std::size_t spin, std::size_t kpoint, Kpoint_data<T>& data) const
{
	assert(spin < n_spins_);
	assert(kpoint < n_kpoints_);
	assert(is_single_precision() == (std::is_same<T, float>::value));
	
	data.energies.resize(n_bands_);
	data.occupations.resize(n_bands_);

	auto record = 2 + (n_bands_ + 1) * (spin * n_kpoints_ + kpoint);
	seek_record(record);

	double n_plane_waves;
	read(n_plane_waves);
	data.n_plane_waves = to_positive_sizet(n_plane_waves);

	read(data.k);

	for (std::size_t i = 0; i < n_bands_; ++i)
	{
		read(data.energies[i]);
		file_.ignore(8);				// Skip the energy imaginary part, should be zero
		read(data.occupations[i]);
	}

	data.gs.reserve(data.n_plane_waves);
	compute_g_lattice(data.k, data.gs);

	if (data.gs.size() != data.n_plane_waves)
		throw Bad_wavecar_file("Inconsistent number of plane waves");

	data.coeffs.resize(data.n_plane_waves, n_bands_);
	for (std::size_t i = 0; i < n_bands_; ++i)
	{
		seek_record(++record);
		read(&data.coeffs(0, i), data.n_plane_waves);
	}
}

template void Wavecar_reader::get_kpoint_data(
	std::size_t, std::size_t, Kpoint_data<float>&) const;

template void Wavecar_reader::get_kpoint_data(
	std::size_t, std::size_t, Kpoint_data<double>&) const;

//////////////////////////////////////////////////////////////////////////

void Wavecar_reader::read_header()
{
	double record_length, n_spins, r_tag, n_kpoints, n_bands;

	read(record_length);
	read(n_spins);
	read(r_tag);

	record_length_ = to_positive_sizet(record_length);
	n_spins_ = to_positive_sizet(n_spins);

	if (r_tag == single_precision_tag)
		precision_ = Precision::SINGLE;
	else if (r_tag == double_precision_tag)
		precision_ = Precision::DOUBLE;
	else
		throw Bad_wavecar_file("Unsupported RTAG value");

	seek_record(1);

	read(n_kpoints);
	n_kpoints_ = to_positive_sizet(n_kpoints);

	read(n_bands);
	n_bands_ = to_positive_sizet(n_bands);

	read(e_cut_);
	read(a0_);
	read(a1_);
	read(a2_);
}

void Wavecar_reader::compute_reciprocal()
{
	const auto uc_volume = det(a0_, a1_, a2_);
	b0_ = 2 * PI / uc_volume * cross(a1_, a2_);
	b1_ = 2 * PI / uc_volume * cross(a2_, a0_);
	b2_ = 2 * PI / uc_volume * cross(a0_, a1_);

	const double g_max_over_2pi = std::sqrt(TWO_M_OVER_HBAR_SQ * e_cut_) / (2 * PI);

	// NB: for oblique unit cell (i_m) is not (Gm / |b_i|), but (Gm * |a_i| / 2pi)
	max_g0_ = static_cast<std::size_t>(std::floor(g_max_over_2pi * a0_norm())) + 1;
	max_g1_ = static_cast<std::size_t>(std::floor(g_max_over_2pi * a1_norm())) + 1;
	max_g2_ = static_cast<std::size_t>(std::floor(g_max_over_2pi * a2_norm())) + 1;
}

void Wavecar_reader::compute_g_lattice(
	const Vec3<double>& k, std::vector<Vec3<std::size_t>>& gs) const
{
	const auto two_m_e_cut_over_hbar_sq = TWO_M_OVER_HBAR_SQ * e_cut_;

	gs.clear();
	for (std::size_t i2 = 0; i2 < size_g2(); ++i2)
	{
		const auto i2s = index_shift(i2, max_g2_);
		const auto g2 = (k[2] + i2s) * b2_;
		for (std::size_t i1 = 0; i1 < size_g1(); ++i1)
		{
			const auto i1s = index_shift(i1, max_g1_);
			const auto g2_p_g1 = g2 + (k[1] + i1s) * b1_;
			for (std::size_t i0 = 0; i0 < size_g0(); ++i0)
			{
				const auto i0s = index_shift(i0, max_g0_);
				const auto g = g2_p_g1 + (k[0] + i0s) * b0_;
				const auto norm_g_sq = norm_sq(g);
				if (norm_g_sq < two_m_e_cut_over_hbar_sq)
					gs.push_back({i0, i1, i2});
			}
		}
	}
}

void Wavecar_reader::seek_record(std::size_t n) const
{
	file_.seekg(static_cast<unsigned long long>(n) * record_length_);
}

std::size_t Wavecar_reader::to_positive_sizet(double x)
{
	auto r = static_cast<std::size_t>(x);
	if (x <= 0 || x != r)
		throw Bad_wavecar_file("Positive integral value expected");

	return r;
}
