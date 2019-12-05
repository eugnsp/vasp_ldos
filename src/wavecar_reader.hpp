#pragma once
#include "matrix.hpp"
#include "vec3.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

template<typename T>
struct Kpoint_data
{
	Vec3<double> k;
	std::size_t n_plane_waves;
	std::vector<double> energies;
	std::vector<double> occupations;
	Matrix<std::complex<T>> coeffs;
	std::vector<Vec3<std::size_t>> gs;
};

class Wavecar_reader
{
public:
	Wavecar_reader(const std::string& filename)
	{
		file_.exceptions(std::ifstream::badbit | std::ifstream::failbit);
		file_.open(filename, std::ifstream::binary);

		read_header();
		compute_reciprocal();
	}

	bool is_single_precision() const
	{
		return precision_ == Precision::SINGLE;
	}

	bool is_double_precision() const
	{
		return precision_ == Precision::DOUBLE;
	}

	std::size_t n_spins() const
	{
		return n_spins_;
	}

	std::size_t n_kpoints() const
	{
		return n_kpoints_;
	}

	std::size_t n_bands() const
	{
		return n_bands_;
	}

	double e_cut() const
	{
		return e_cut_;
	}

	const Vec3<double>& a0() const
	{
		return a0_;
	}

	const Vec3<double>& a1() const
	{
		return a1_;
	}

	const Vec3<double>& a2() const
	{
		return a2_;
	}

	double a0_norm() const
	{
		return norm(a0_);
	}

	double a1_norm() const
	{
		return norm(a1_);
	}

	double a2_norm() const
	{
		return norm(a2_);
	}

	std::size_t size_g0() const
	{
		return 2 * max_g0_ + 1;
	}

	std::size_t size_g1() const
	{
		return 2 * max_g1_ + 1;
	}

	std::size_t size_g2() const
	{
		return 2 * max_g2_ + 1;
	}

	template<typename T>
	void get_kpoint_data(std::size_t spin, std::size_t kpoint, Kpoint_data<T>& data)
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
			throw std::runtime_error("Bad WAVECAR: Inconsistent number of plane waves");

		data.coeffs.resize(data.n_plane_waves, n_bands_);
		for (std::size_t i = 0; i < n_bands_; ++i)
		{
			seek_record(++record);
			read(&data.coeffs(0, i), data.n_plane_waves);
		}
	}

private:
	enum class Precision
	{
		SINGLE,
		DOUBLE
	};

	static constexpr double PI = 3.141592653589793238463;
	static constexpr double TWO_M_OVER_HBAR_SQ = 0.262465831;
	static constexpr double single_precision_tag = 45'200;
	static constexpr double double_precision_tag = 45'210;

	template<typename T>
	void read(T& x)
	{
		file_.read(reinterpret_cast<char*>(&x), sizeof(T));
	}

	template<typename T>
	void read(T* buff, std::size_t count)
	{
		file_.read(reinterpret_cast<char*>(buff), sizeof(T) * count);
	}

	void read_header()
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
			throw std::runtime_error("Bad WAVECAR: Unsupported RTAG value");

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

	void compute_reciprocal()
	{
		const auto uc_volume = a0_ * (a1_ ^ a2_);
		b0_ = 2 * PI / uc_volume * (a1_ ^ a2_);
		b1_ = 2 * PI / uc_volume * (a2_ ^ a0_);
		b2_ = 2 * PI / uc_volume * (a0_ ^ a1_);

		const double g_max_over_2pi = std::sqrt(TWO_M_OVER_HBAR_SQ * e_cut_) / (2 * PI);

		// NB: for oblique unit cell (i_m) is not (Gm / |b_i|), but (Gm * |a_i| / 2pi)
		max_g0_ = static_cast<std::size_t>(std::floor(g_max_over_2pi * a0_norm())) + 1;
		max_g1_ = static_cast<std::size_t>(std::floor(g_max_over_2pi * a1_norm())) + 1;
		max_g2_ = static_cast<std::size_t>(std::floor(g_max_over_2pi * a2_norm())) + 1;
	}

	void compute_g_lattice(const Vec3<double>& k, std::vector<Vec3<std::size_t>>& gs) const
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

	void seek_record(std::size_t n)
	{
		file_.seekg(static_cast<unsigned long long>(n) * record_length_);
	}

	// Casts a double to a (positive) integer checking whether
	// it can be represented as a positive integer
	static std::size_t to_positive_sizet(double x)
	{
		const auto r = static_cast<std::size_t>(x);
		if (x <= 0 || x != r)
			throw std::runtime_error("Bad WAVECAR: Positive integral value expected");

		return r;
	}

	template<typename T>
	static typename std::make_signed_t<T> index_shift(T i, T i_max)
	{
		using Ts = std::make_signed_t<T>;

		const auto i_max_s = static_cast<Ts>(i_max);
		auto index = static_cast<Ts>(i);
		if (index > i_max_s)
			index -= (2 * i_max_s + 1);

		return index;
	}

private:
	std::size_t record_length_;

	std::size_t n_spins_;
	std::size_t n_kpoints_;
	std::size_t n_bands_;

	// Maximum index (i_m) of reciprocal lattice vectors in the expansion,
	// the total number of vectors is (2 * i_m + 1)
	std::size_t max_g0_;
	std::size_t max_g1_;
	std::size_t max_g2_;

	double e_cut_;
	Vec3<double> a0_, a1_, a2_;		// Direct lattice
	Vec3<double> b0_, b1_, b2_;		// Reciprocal lattice

	Precision precision_;

	std::ifstream file_;
};
