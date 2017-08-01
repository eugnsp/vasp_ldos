#pragma once
#include "matrix.hpp"
#include <cstddef>
#include <complex>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

template<typename T>
using Vec3 = std::array<T, 3>;

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
	explicit Wavecar_reader(const std::string& filename);

	bool is_single_precision() const;
	bool is_double_precision() const;

	std::size_t n_spins() const;
	std::size_t n_kpoints() const;
	std::size_t n_bands() const;
	double e_cut() const;

	const Vec3<double>& a0() const;
	const Vec3<double>& a1() const;
	const Vec3<double>& a2() const;

	double a0_norm() const;
	double a1_norm() const;
	double a2_norm() const;

	std::size_t size_g0() const;
	std::size_t size_g1() const;
	std::size_t size_g2() const;

	template<typename T>
	void get_kpoint_data(std::size_t spin, std::size_t kpoint, Kpoint_data<T>&) const;

private:
	void read_header();
	void compute_reciprocal();

	void compute_g_lattice(
		const Vec3<double>& k, std::vector<Vec3<std::size_t>>&) const;

	void seek_record(std::size_t) const;

	// Skips sizeof(T) bytes in the file stream
	template<typename T>
	void skip() const;

	template<typename T>
	void read(T&) const;

	template<typename T>
	void read(T*, std::size_t n_elements) const;

	// Casts a double to a (positive) integer checking whether
	// it can be represented as a positive integer
	static std::size_t to_positive_sizet(double);

private:
	enum class Precision
	{
		SINGLE,
		DOUBLE
	};

	static constexpr double single_precision_tag = 45200;
	static constexpr double double_precision_tag = 45210;

private:
	mutable std::ifstream file_;
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
};

class Bad_wavecar_file : public std::runtime_error
{
public:
	Bad_wavecar_file(const std::string& err)
		: std::runtime_error("Bad WAVECAR file: " + err)
	{ }
};
