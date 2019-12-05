#pragma once
#include <algorithm>
#include <vector>
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
		: rows_(rows), cols_(cols)
	{
		data_.resize(rows_ * cols_);
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
		data_.resize(rows_ * cols_);
	}

	T& operator()(std::size_t row, std::size_t col)
	{
		assert(row < rows_);
		assert(col < cols_);

		return data_[row + col * rows_];
	}

	const T& operator()(std::size_t row, std::size_t col) const
	{
		assert(row < rows_);
		assert(col < cols_);

		return data_[row + col * rows_];
	}

	T* data()
	{
		return data_.data();
	}

	const T* data() const
	{
		return data_.data();
	}

	void fill(const T& value)
	{
		std::fill(data_.begin(), data_.end(), value);
	}

private:
	std::vector<T> data_;

	std::size_t rows_ = 0;
	std::size_t cols_ = 0;
};
