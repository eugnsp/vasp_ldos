#pragma once
#include <fftw3.h>

#include <cassert>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

template<typename T>
class Fft
{
	static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, "Bad data type");

public:
	Fft(std::size_t size, std::size_t n_transforms, std::complex<T>* data)
	{
		assert(size > 0);
		assert(n_transforms > 0);
		const int n = static_cast<int>(size);

		if constexpr (std::is_same_v<T, float>)
			plan_ = fftwf_plan_many_dft(1, &n, static_cast<int>(n_transforms),
				reinterpret_cast<fftwf_complex*>(data), nullptr, 1, n,
				reinterpret_cast<fftwf_complex*>(data), nullptr, 1, n, FFTW_BACKWARD, FFTW_ESTIMATE);
		else
			plan_ = fftw_plan_many_dft(1, &n, static_cast<int>(n_transforms),
				reinterpret_cast<fftw_complex*>(data), nullptr, 1, n,
				reinterpret_cast<fftw_complex*>(data), nullptr, 1, n, FFTW_BACKWARD, FFTW_ESTIMATE);

		if (!plan_)
			throw std::runtime_error("FFTW plan creation failed");
	}

	~Fft()
	{
		if (plan_)
		{
			if constexpr (std::is_same_v<T, float>)
			{
				fftwf_destroy_plan(plan_);
				fftwf_cleanup();
			}
			else
			{
				fftw_destroy_plan(plan_);
				fftw_cleanup();
			}
		}
	}

	Fft& operator=(const Fft&) = delete;

	void transform() const
	{
		if constexpr (std::is_same_v<T, float>)
			fftwf_execute(plan_);
		else
			fftw_execute(plan_);
	}

private:
	std::conditional_t<std::is_same_v<T, float>, fftwf_plan, fftw_plan> plan_ = nullptr;
};
