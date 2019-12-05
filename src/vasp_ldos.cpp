#include "command_line.hpp"
#include "ldos_writer.hpp"
#include "matrix.hpp"
#include "wavecar_reader.hpp"

#ifdef USE_MKL_FFT
	#include "fft_mkl.hpp"
#else
	#include "fft_fftw.hpp"
#endif

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

enum class Cell_direction
{
	A0, A1, A2
};

struct Fft_size
{
	std::size_t size;
	std::size_t n_transforms;
};

Cell_direction get_direction(const Wavecar_reader& reader)
{
	if (reader.a0_norm() > reader.a1_norm() && reader.a0_norm() > reader.a2_norm())
		return Cell_direction::A0;
	else if (reader.a1_norm() > reader.a2_norm() && reader.a1_norm() > reader.a0_norm())
		return Cell_direction::A1;
	else if (reader.a2_norm() > reader.a0_norm() && reader.a2_norm() > reader.a1_norm())
		return Cell_direction::A2;
	else
		throw std::runtime_error("Bad supercell size");
}

double get_height(const Wavecar_reader& wc_reader, Cell_direction dir)
{
	switch (dir)
	{
	case Cell_direction::A0:
		return wc_reader.a0_norm();

	case Cell_direction::A1:
		return wc_reader.a1_norm();

	default: // case Cell_direction::A2:
		return wc_reader.a2_norm();
	}
}

Fft_size get_fft_size(Wavecar_reader& wc_reader, Cell_direction dir)
{
	switch (dir)
	{
	case Cell_direction::A0:
		return {wc_reader.size_g0(), wc_reader.size_g1() * wc_reader.size_g2()};

	case Cell_direction::A1:
		return {wc_reader.size_g1(), wc_reader.size_g2() * wc_reader.size_g0()};

	default: // case Cell_direction::A2:
		return {wc_reader.size_g2(), wc_reader.size_g0() * wc_reader.size_g1()};
	}
}

template<typename T>
void map_g_sphere_to_fft_blocks(const Wavecar_reader& wc_reader,
								Matrix<std::complex<T>>& cs,
                                const Kpoint_data<T>& kpoint_data,
								std::size_t band, Cell_direction dir)
{
	cs.fill(0);

	switch (dir)
	{
	case Cell_direction::A0:
		for (std::size_t ipw = 0; ipw < kpoint_data.n_plane_waves; ++ipw)
		{
			const auto& g = kpoint_data.gs[ipw];
			const auto g_parallel_index = g[1] + g[2] * wc_reader.size_g1();
			cs(g[0], g_parallel_index) = kpoint_data.coeffs(ipw, band);
		}
		break;

	case Cell_direction::A1:
		for (std::size_t ipw = 0; ipw < kpoint_data.n_plane_waves; ++ipw)
		{
			const auto& g = kpoint_data.gs[ipw];
			const auto g_parallel_index = g[2] + g[0] * wc_reader.size_g2();
			cs(g[1], g_parallel_index) = kpoint_data.coeffs(ipw, band);
		}
		break;

	case Cell_direction::A2:
		for (std::size_t ipw = 0; ipw < kpoint_data.n_plane_waves; ++ipw)
		{
			const auto& g = kpoint_data.gs[ipw];
			const auto g_parallel_index = g[0] + g[1] * wc_reader.size_g0();
			cs(g[2], g_parallel_index) = kpoint_data.coeffs(ipw, band);
		}
	}
}

template<typename T>
void process(Wavecar_reader& reader, Ldos_writer& writer, Cell_direction dir)
{
	const auto fft_size = get_fft_size(reader, dir);

	Matrix<std::complex<T>> cs{fft_size.size, fft_size.n_transforms};
	Matrix<float> cs_sq{fft_size.size, reader.n_bands()};
	Kpoint_data<T> kpoint_data;

	Fft<T> fft{fft_size.size, fft_size.n_transforms, cs.data()};

	auto energy_min = std::numeric_limits<double>::max();
	auto energy_max = -std::numeric_limits<double>::max();
	auto cs_sq_max = -std::numeric_limits<float>::max();

	std::cout << std::string(reader.n_spins() * reader.n_kpoints(), '*') << std::endl;

	for (std::size_t is = 0; is < reader.n_spins(); ++is)
		for (std::size_t ik = 0; ik < reader.n_kpoints(); ++ik)
		{
			reader.get_kpoint_data(is, ik, kpoint_data);

			cs_sq.fill(0);
			for (std::size_t ib = 0; ib < reader.n_bands(); ++ib)
			{
				energy_min = std::min(energy_min, kpoint_data.energies[ib]);
				energy_max = std::max(energy_max, kpoint_data.energies[ib]);

				map_g_sphere_to_fft_blocks(reader, cs, kpoint_data, ib, dir);
				fft.transform();

				// Sum over G||
				for (std::size_t ip = 0; ip < fft_size.n_transforms; ++ip)
					for (std::size_t il = 0; il < fft_size.size; ++il)
					{
						const auto sq = static_cast<float>(std::norm(cs(il, ip)));
						cs_sq_max = std::max(cs_sq_max, sq);
						cs_sq(il, ib) += sq;
					}
			}

			writer.write_ldos(kpoint_data.k, kpoint_data.energies, kpoint_data.occupations, cs_sq);
			std::cout << '.' << std::flush;
		}

	writer.write_minmax_values(energy_min, energy_max, cs_sq_max);
	std::cout << std::endl;
}

void print_wavecar_info(const Wavecar_reader& reader)
{
	std::cout << "WAVECAR file:\n"
			  << "Precision: " << (reader.is_single_precision() ? "single" : "double") << '\n'
			  << "Number of spin components: " << reader.n_spins() << '\n'
			  << "Number of k-points: " << reader.n_kpoints() << '\n'
			  << "Number of bands: " << reader.n_bands() << '\n'
			  << "Cut-off energy: " << reader.e_cut() << " eV\n\n"

			  << std::fixed << std::setprecision(5) << "Direct lattice:\n"
			  << " a1 = " << reader.a()[0] << " Ang\n"
			  << " a2 = " << reader.a()[1] << " Ang\n"
			  << " a3 = " << reader.a()[2] << " Ang\n\n"

			  << "G-lattice size: " << reader.size_g0() << " x "
			  << reader.size_g1() << " x " << reader.size_g2() << '\n' << std::endl;
}

void print_help()
{
	std::cout << "Synopsis:\n"
			  << "    vasp_ldos [options]\n"
			  << "Options:\n"
			  << "    -h               print help\n"
			  << "    -o <name>        output LDOS filename (no default)\n"
			  << "    -w <name>        input WAVECAR filename (default: \"WAVECAR\")\n"
			  << "    -f <value>       Fermi level value (default: 0)\n"
			  << "    -c <comment>     arbitrary text comment (default: none)\n\n"
			  << "If no output filename is given, WAVECAR file basic\n"
			  << "information is displayed and the program terminates." << std::endl;
}

//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	try
	{
		const Command_line cl(argc, argv);

		if (cl.option_exists("-h"))
		{
			print_help();
			return 0;
		}

		const std::string wc_filename = cl.get_option_or("-w", "WAVECAR");
		Wavecar_reader reader(wc_filename);
		print_wavecar_info(reader);

		if (!cl.option_exists("-o"))
			return 0;

		const std::string output_filename = cl.get_option("-o");
		const auto user_comment = cl.get_option_or("-c", "");
		const double fermi_energy = std::stod(cl.get_option_or("-f", "0"));

		const auto cell_direction = get_direction(reader);
		Ldos_writer writer(output_filename, reader, get_fft_size(reader, cell_direction).size,
			get_height(reader, cell_direction), fermi_energy, user_comment);

		if (reader.is_single_precision())
			process<float>(reader, writer, cell_direction);
		else
			process<double>(reader, writer, cell_direction);
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception!\n" << e.what() << std::endl;
		return -1;
	}
	catch (...)
	{
		std::cerr << "Exception!" << std::endl;
		return -1;
	}

	std::cout << "Done!" << std::endl;
	return 0;
}
