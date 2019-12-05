#include "ldos_writer.hpp"
#include <cstdint>
#include <string>
#include <ctime>
#include <cassert>

namespace
{
std::string date_time_string(const char* format = "%a, %d %b %Y %T")
{
	std::time_t tm = std::time(nullptr);
	std::tm ltm;
#ifdef __unix__
	localtime_r(&tm, &ltm);
#else
	localtime_s(&ltm, &tm);
#endif

	char time[100];
	std::strftime(time, 100, format, &ltm);
	return std::string(time);
}
}

Ldos_writer::Ldos_writer(
	const std::string& filename,
	std::size_t n_spins, std::size_t n_kpoints,
	std::size_t n_bands, std::size_t n_layers,
	double supercell_height, double fermi_energy,
	const std::string& user_comment /* = {} */)
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

void Ldos_writer::write_ldos(
	const Vec3& k, const std::vector<double>& energies,
	const std::vector<double>& occupations, const Matrix<float>& cs_sq)
{
	assert(energies.size() == n_bands_ && occupations.size() == n_bands_);
	assert(cs_sq.rows() == n_layers_ && cs_sq.cols() == n_bands_);

	write(k);
	write(energies.data(), energies.size());
	write(occupations.data(), occupations.size());
	write(cs_sq.data(), cs_sq.size());
}

void Ldos_writer::write_minmax_values(
	double energy_min, double energy_max, float cs_sq_max)
{
	assert(energy_min < energy_max);

	file_.seekp(minmax_values_pos_);
 	write(energy_min);
 	write(energy_max);
 	write(cs_sq_max);
}

template<typename T>
void Ldos_writer::write(const T& x)
{
	file_.write(reinterpret_cast<const char*>(&x), sizeof(T));
}

template<typename T>
void Ldos_writer::write(const T* x, std::size_t n_elements)
{
	file_.write(reinterpret_cast<const char*>(x), sizeof(T) * n_elements);
}
