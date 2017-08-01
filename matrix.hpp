#pragma once
#include <malloc.h>
#include <cstddef>
#include <cstring>
#include <new>
#include <cassert>

// Simple class for column-major matricies
template<typename T>
class Matrix
{
public:
	Matrix() = default;

	Matrix(std::size_t rows, std::size_t cols)
		: rows_(rows), cols_(cols), capacity_(rows * cols)
	{
		assert(rows > 0 && cols > 0);
		allocate();
	}

	~Matrix()
	{
		deallocate();
	}

	std::size_t rows() const
	{
		return rows_;
	}

	std::size_t cols() const
	{
		return cols_;
	}

	std::size_t size() const
	{
		return rows_ * cols_;
	}

	void resize(std::size_t rows, std::size_t cols)
	{
		assert(rows > 0 && cols > 0);

		rows_ = rows;
		cols_ = cols;
		if (capacity_ < rows_ * cols_)
		{
			deallocate();
			allocate();
			capacity_ = rows_ * cols_;
		}
	}

	T& operator()(std::size_t row, std::size_t col)
	{
		assert(row < rows_);
		assert(col < cols_);
		
		return data_[row + col * rows_];
	}

	T operator()(std::size_t row, std::size_t col) const
	{
		assert(row < rows_);
		assert(col < cols_);

		return data_[row + col * rows_];
	}

	T* data()
	{
		return data_;
	}

	const T* data() const
	{
		return data_;
	}

	void zero()
	{
		std::memset(data_, 0, size() * sizeof(T));
	}

private:
	void allocate()
	{
		if (size() > static_cast<std::size_t>(-1) / sizeof(T))
			throw std::bad_array_new_length();

		data_ = reinterpret_cast<T*>(_mm_malloc(size() * sizeof(T), 64));
		if (!data_)
			throw std::bad_alloc();
	}

	void deallocate()
	{
		_mm_free(data_);
	}

private:
	T* data_ = nullptr;

	std::size_t rows_ = 0;
	std::size_t cols_ = 0;
	std::size_t capacity_ = 0;
};
