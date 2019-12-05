#pragma once
#include "matrix.hpp"
#include "vec3.hpp"

#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

class Ldos_writer
{
public:
	Ldos_writer(const std::string& filename,
				std::size_t n_spins, std::size_t n_kpoints,
				std::size_t n_bands, std::size_t n_layers,
				double supercell_height, double fermi_energy,
				const std::string& user_comment)
		: n_bands_(n_bands), n_layers_(n_layers)
	{
		assert(n_spins > 0);
		assert(n_kpoints > 0);
		assert(n_bands > 0);
		assert(n_layers > 0);

		file_.exceptions(std::ofstream::failbit | std::ofstream::badbit);
		file_.open(filename, std::ofstream::binary);

		const std::size_t header_length = 500;
		std::string header("Depth-k resolved DOS data file, created on: ");
		header += date_time_string() + "; ";
		header += std::to_string(n_kpoints) + " k points, " +
			std::to_string(n_bands) + " bands, " + std::to_string(n_layers) + " layers";

		if (!user_comment.empty())
			header += "; Comment: " + user_comment;

		header.resize(header_length, ' ');
		write(header.c_str(), header.length());

		const std::uint32_t file_format_version = 103;
		write(file_format_version);

		write(static_cast<std::uint32_t>(n_spins));
		write(static_cast<std::uint32_t>(n_kpoints));
		write(static_cast<std::uint32_t>(n_bands));
		write(static_cast<std::uint32_t>(n_layers));
		write(supercell_height);
		write(fermi_energy);

		minmax_values_pos_ = file_.tellp();
		write(0.);	// Reserved for energy_min
		write(0.);	// Reserved for energy_max
		write(0.f); // Reserved for cs_sq_max
	}

	void write_ldos(const Vec3<double>& k, const std::vector<double>& energies,
		  			const std::vector<double>& occupations,
					const Matrix<float>& cs_sq)
	{
		assert(energies.size() == n_bands_ && occupations.size() == n_bands_);
		assert(cs_sq.rows() == n_layers_ && cs_sq.cols() == n_bands_);

		write(k);
		write(energies.data(), energies.size());
		write(occupations.data(), occupations.size());
		write(cs_sq.data(), cs_sq.size());
	}

	void write_minmax_values(double energy_min, double energy_max, float cs_sq_max)
	{
		assert(energy_min < energy_max);

		file_.seekp(minmax_values_pos_);
	 	write(energy_min);
	 	write(energy_max);
	 	write(cs_sq_max);
	}

private:
	template<typename T>
	void write(const T& x)
	{
		file_.write(reinterpret_cast<const char*>(&x), sizeof(T));
	}

	template<typename T>
	void write(const T* buff, std::size_t count )
	{
		file_.write(reinterpret_cast<const char*>(buff), sizeof(T) * count);
	}

	static std::string date_time_string()
	{
		const auto now = std::chrono::system_clock::now();
    	const auto time = std::chrono::system_clock::to_time_t(now);

    	std::stringstream ss;
    	ss << std::put_time(std::localtime(&time), "%a, %d %b %Y %T");
		return ss.str();
	}

private:
	std::ofstream file_;
	std::streampos minmax_values_pos_;

	const std::size_t n_bands_;
	const std::size_t n_layers_;
};
