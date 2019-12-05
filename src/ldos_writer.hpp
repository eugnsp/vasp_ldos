#pragma once
#include "matrix.hpp"
#include <cstddef>
#include <array>
#include <vector>
#include <string>
#include <fstream>

class Ldos_writer
{
private:
	using Vec3 = std::array<double, 3>;

public:
	Ldos_writer(
		const std::string& filename,		
		std::size_t n_spins, std::size_t n_kpoints,
		std::size_t n_bands, std::size_t n_layers,
		double supercell_height, double e_fermi,
		const std::string& user_comment = {});

	void write_ldos(
		const Vec3& k, const std::vector<double>& energies,
		const std::vector<double>& occupations,
		const Matrix<float>& cs_sq);

	void write_minmax_values(double energy_min, double energy_max, float cs_sq_max);

private:
	template<typename T>
	void write(const T&);

	template<typename T>
	void write(const T*, std::size_t n_elements);

private:
	std::ofstream file_;
	std::streampos minmax_values_pos_;

	const std::size_t n_bands_;
	const std::size_t n_layers_;
};
