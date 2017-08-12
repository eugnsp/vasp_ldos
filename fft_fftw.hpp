#pragma once
#include <fftw3.h>
#include <cstddef>
#include <complex>
#include <type_traits>
#include <stdexcept>
#include <cassert>

template<typename T>
class Fft
{
	static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value,
				  "Bad data type");

public:
	Fft(std::size_t size, std::size_t n_transforms, std::complex<T>*);
	~Fft();

	void transform() const;

private:
	std::conditional_t<std::is_same<T, float>::value,
		fftwf_plan, fftw_plan> plan_ = nullptr;
};

//////////////////////////////////////////////////////////////////////////

inline Fft<float>::Fft(std::size_t size, std::size_t n_transforms, std::complex<float>* data)
{
	assert(size > 0);
	assert(n_transforms > 0);

	const int n = static_cast<int>(size);
	plan_ = fftwf_plan_many_dft(
		1, &n, static_cast<int>(n_transforms),
		reinterpret_cast<fftwf_complex*>(data), nullptr, 1, n, 
		reinterpret_cast<fftwf_complex*>(data),	nullptr, 1, n,
		FFTW_BACKWARD, FFTW_ESTIMATE);

	if (!plan_)
		throw std::runtime_error("FFTW plan creation failed");
}

inline Fft<double>::Fft(std::size_t size, std::size_t n_transforms, std::complex<double>* data)
{
	assert(size > 0);
	assert(n_transforms > 0);

	const int n = static_cast<int>(size);
	plan_ = fftw_plan_many_dft(
		1, &n, static_cast<int>(n_transforms),
		reinterpret_cast<fftw_complex*>(data), nullptr, 1, n,
		reinterpret_cast<fftw_complex*>(data), nullptr, 1, n,
		FFTW_BACKWARD, FFTW_ESTIMATE);

	if (!plan_)
		throw std::runtime_error("FFTW plan creation failed");
}

inline Fft<float>::~Fft()
{
	if (plan_)
	{
		fftwf_destroy_plan(plan_);
		fftwf_cleanup();
	}
}

inline Fft<double>::~Fft()
{
	if (plan_)
	{
		fftw_destroy_plan(plan_);
		fftw_cleanup();
	}
}

inline void Fft<float>::transform() const
{
	fftwf_execute(plan_);
}

inline void Fft<double>::transform() const
{
	fftw_execute(plan_);
}
