#include "command_line.hpp"
#include "matrix.hpp"
#include "fft_mkl.hpp"
//#include "fft_fftw.hpp"
#include "wavecar_reader.hpp"
#include "outcar_reader.hpp"
#include "ldos_writer.hpp"
#include <cstddef>
#include <utility>
#include <limits>
#include <complex>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <memory>
#include <type_traits>
#include <exception>
#include <stdexcept>
#include <cassert>

template<typename T, typename... Args>
std::unique_ptr<T> my_make_unique(Args&&... args)
{
	return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

class App
{
public:
	App(int argc, char* argv[]);

	void run();

private:
	enum class Height_direction
	{
		A0, A1, A2
	};
	
	struct Fft_size
	{
		std::size_t size;
		std::size_t n_transforms;
	};

	template<typename T>
	void process();

	double get_height() const;
	Fft_size get_fft_size() const;

	template<typename T>
	void map_g_sphere_to_fft_blocks(
		Matrix<std::complex<T>>&, const Kpoint_data<T>&, std::size_t band) const;

	void print_outcar_info() const;
	void print_wavecar_info() const;
	static void print_help();

private:
	std::unique_ptr<Wavecar_reader> wavecar_reader_;
	std::unique_ptr<Outcar_reader> outcar_reader_;
	std::unique_ptr<Ldos_writer> writer_;
	Height_direction height_direction_;

	const Command_line cl;
};

App::App(int argc, char* argv[])
	: cl(argc, argv)
{ }

void App::run()
{
	if (cl.option_exists("-h"))
	{
		print_help();
		return;
	}

	const std::string outcar_filename = cl.get_option("-o", "OUTCAR");
	outcar_reader_ = my_make_unique<Outcar_reader>(outcar_filename);
	print_outcar_info();

	const std::string wavecar_filename = cl.get_option("-w", "WAVECAR");
	wavecar_reader_ = my_make_unique<Wavecar_reader>(wavecar_filename);
	print_wavecar_info();

	if (!cl.option_exists("-l"))
		return;

	if (wavecar_reader_->a0_norm() > wavecar_reader_->a1_norm() && wavecar_reader_->a0_norm() > wavecar_reader_->a2_norm())
		height_direction_ = Height_direction::A0;
	else if (wavecar_reader_->a1_norm() > wavecar_reader_->a2_norm() && wavecar_reader_->a1_norm() > wavecar_reader_->a0_norm())
		height_direction_ = Height_direction::A1;
	else if (wavecar_reader_->a2_norm() > wavecar_reader_->a0_norm() && wavecar_reader_->a2_norm() > wavecar_reader_->a1_norm())
		height_direction_ = Height_direction::A2;
	else
		throw std::runtime_error("Bad supercell size");

	const std::string ldos_filename = cl.get_option("-l");
	std::string user_comment;
	if (cl.option_exists("-c"))
		user_comment = cl.get_option("-c");

	writer_ = my_make_unique<Ldos_writer>(
		ldos_filename, wavecar_reader_->n_spins(), wavecar_reader_->n_kpoints(),
		wavecar_reader_->n_bands(), get_fft_size().size, get_height(), outcar_reader_->fermi_energy(),
		user_comment);

	if (wavecar_reader_->is_single_precision())
		process<float>();
	else
		process<double>();
}

template<typename T>
void App::process()
{
	const auto fft_size = get_fft_size();

	Matrix<std::complex<T>> cs{fft_size.size, fft_size.n_transforms};
	Matrix<float> cs_sq{fft_size.size, wavecar_reader_->n_bands()};
	Kpoint_data<T> kpoint_data;

	Fft<T> fft{fft_size.size, fft_size.n_transforms, cs.data()};

	auto energy_min = std::numeric_limits<double>::max();
	auto energy_max = -std::numeric_limits<double>::max();
	auto cs_sq_max = -std::numeric_limits<float>::max();

	std::cout << std::string(wavecar_reader_->n_spins() * wavecar_reader_->n_kpoints(), '*') << '\n' << std::flush;

	for (std::size_t is = 0; is < wavecar_reader_->n_spins(); ++is)
		for (std::size_t ik = 0; ik < wavecar_reader_->n_kpoints(); ++ik)
		{
			std::cout << '.' << std::flush;
			wavecar_reader_->get_kpoint_data(is, ik, kpoint_data);

			const auto n_plane_waves = kpoint_data.gs.size();
			cs_sq.zero();

			for (std::size_t ib = 0; ib < wavecar_reader_->n_bands(); ++ib)
			{
				energy_min = std::min(energy_min, kpoint_data.energies[ib]);
				energy_max = std::max(energy_max, kpoint_data.energies[ib]);

				map_g_sphere_to_fft_blocks(cs, kpoint_data, ib);
				fft.transform();

				// Sum over G||
				for (std::size_t ip = 0; ip < fft_size.n_transforms; ++ip)
					for (std::size_t il = 0; il < fft_size.size; ++il)
					{
						const auto c = cs(il, ip);

						const auto re = c.real();
						const auto im = c.imag();
						const auto sq = static_cast<float>(re * re + im * im);
						cs_sq_max = std::max(cs_sq_max, sq);
						cs_sq(il, ib) += sq;
					}
			}

			writer_->write_ldos(
				kpoint_data.k, kpoint_data.energies, kpoint_data.occupations, cs_sq);
		}

	std::cout << "\n\n";
	writer_->write_minmax_values(energy_min, energy_max, cs_sq_max);
}

double App::get_height() const
{
	switch (height_direction_)
	{
	case Height_direction::A0:
		return wavecar_reader_->a0_norm();

	case Height_direction::A1:
		return wavecar_reader_->a1_norm();

	default: // case Z_direction::A2:
		return wavecar_reader_->a2_norm();
	}
}

auto App::get_fft_size() const -> Fft_size
{
	switch (height_direction_)
	{
	case Height_direction::A0:
		return {wavecar_reader_->size_g0(), wavecar_reader_->size_g1() * wavecar_reader_->size_g2()};

	case Height_direction::A1:
		return {wavecar_reader_->size_g1(), wavecar_reader_->size_g2() * wavecar_reader_->size_g0()};

	default: // A2
		return {wavecar_reader_->size_g2(), wavecar_reader_->size_g0() * wavecar_reader_->size_g1()};
	}
}

template<typename T>
void App::map_g_sphere_to_fft_blocks(
	Matrix<std::complex<T>>& cs, const Kpoint_data<T>& kpoint_data, std::size_t band) const
{
	cs.zero();

	if (height_direction_ == Height_direction::A0)
 		for (std::size_t ipw = 0; ipw < kpoint_data.n_plane_waves; ++ipw)
 		{
 			const auto& g = kpoint_data.gs[ipw];
 			const auto g_parallel_index = g[1] + g[2] * wavecar_reader_->size_g1();
 			cs(g[0], g_parallel_index) = kpoint_data.coeffs(ipw, band);
 		}
 	else if (height_direction_ == Height_direction::A1)
 		for (std::size_t ipw = 0; ipw < kpoint_data.n_plane_waves; ++ipw)
 		{
 			const auto& g = kpoint_data.gs[ipw];
 			const auto g_parallel_index = g[2] + g[0] * wavecar_reader_->size_g2();
 			cs(g[1], g_parallel_index) = kpoint_data.coeffs(ipw, band);
 		}
 	else // if (height_direction_ == Height_direction::A2)
 		for (std::size_t ipw = 0; ipw < kpoint_data.n_plane_waves; ++ipw)
 		{
 			const auto& g = kpoint_data.gs[ipw];
 			const auto g_parallel_index = g[0] + g[1] * wavecar_reader_->size_g0();
 			cs(g[2], g_parallel_index) = kpoint_data.coeffs(ipw, band);
 		}
}

void App::print_outcar_info() const
{ 
	std::cout << "OUTCAR file:\n" << "VASP info: " << outcar_reader_->vasp_info() << '\n'
		<< "Fermi energy: " << outcar_reader_->fermi_energy() << " eV\n\n" << std::flush;
}

void App::print_wavecar_info() const
{
	std::cout << "WAVECAR file:\n"
		<< "Precision: " << (wavecar_reader_->is_single_precision() ? "single" : "double") << '\n'
		<< "Number of spin components: " << wavecar_reader_->n_spins() << '\n'
		<< "Number of k-points: " << wavecar_reader_->n_kpoints() << '\n'
		<< "Number of bands: " << wavecar_reader_->n_bands() << '\n'
		<< "Cut-off energy: " << wavecar_reader_->e_cut() << " eV\n\n"
		
		<< std::fixed << std::setprecision(5)
		<< "Direct lattice:\n"
		<< " a1 = (" << wavecar_reader_->a0()[0] << ", " << wavecar_reader_->a0()[1] << ", " << wavecar_reader_->a0()[2] << ") Ang\n"
		<< " a2 = (" << wavecar_reader_->a1()[0] << ", " << wavecar_reader_->a1()[1] << ", " << wavecar_reader_->a1()[2] << ") Ang\n"
		<< " a3 = (" << wavecar_reader_->a2()[0] << ", " << wavecar_reader_->a2()[1] << ", " << wavecar_reader_->a2()[2] << ") Ang\n\n"

		<< "G-lattice size: "
		<< wavecar_reader_->size_g0() << " x "
		<< wavecar_reader_->size_g1() << " x "
		<< wavecar_reader_->size_g2() << "\n\n" << std::flush;
}

void App::print_help()
{
	std::cout << "Synopsis:\n"
		<< "    ldos [options]\n"
		<< "Options:\n"
		<< "    -h               print help\n"
		<< "    -w <name>        input WAVECAR filename\n"
		<< "    -o <name>        input OUTCAR filename\n"
		<< "    -l <name>        output LDOS filename\n"
		<< "    -c <comment>     arbitrary text comment\n\n"
		<< "If no output filename is given, WAVECAR file basic\n"
		<< "information is displayed and the program terminates.\n";
}

//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	try
	{
		App app(argc, argv);
		app.run();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << "\nError!\n";
		return -1;
	}
	
	std::cout << "Done!\n";
	return 0;
}
