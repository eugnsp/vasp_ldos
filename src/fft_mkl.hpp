#pragma once
#include <mkl_dfti.h>
#include <mkl_service.h>

#include <cassert>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

template<typename T>
class Fft
{
	static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
				  "Bad data type");

public:
	Fft(std::size_t size, std::size_t n_transforms, std::complex<T>* data)
		: data_(data)
	{
		assert(size > 0);
		assert(n_transforms > 0);

		auto status = DftiCreateDescriptor(
			&handle_, (std::is_same_v<T, float> ? DFTI_SINGLE : DFTI_DOUBLE),
			DFTI_COMPLEX, 1, size);

		check_status(status);

		DftiSetValue(handle_, DFTI_COMPLEX_STORAGE, DFTI_COMPLEX_COMPLEX);
		DftiSetValue(handle_, DFTI_PLACEMENT, DFTI_INPLACE);
		DftiSetValue(handle_, DFTI_NUMBER_OF_TRANSFORMS, n_transforms);
		DftiSetValue(handle_, DFTI_INPUT_DISTANCE, size);

		status = DftiCommitDescriptor(handle_);
		check_status(status);
	}

	~Fft()
	{
		if (handle_)
			DftiFreeDescriptor(&handle_);

		mkl_free_buffers();
	}

	Fft& operator=(const Fft&) = delete;

	void transform() const
	{
		auto status = DftiComputeBackward(handle_, data_);
		check_status(status);
	}

private:
	static void check_status(MKL_LONG status)
	{
		if (status && !DftiErrorClass(status, DFTI_NO_ERROR))
			throw std::runtime_error(std::string("MKL FFT failed: ") + DftiErrorMessage(status));
	}

private:
	DFTI_DESCRIPTOR_HANDLE handle_ = nullptr;
	std::complex<T>* const data_;
};
