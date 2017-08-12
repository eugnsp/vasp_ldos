#pragma once
#include <cstddef>
#include <string>
#include <fstream>
#include <stdexcept>

class File_reader
{
protected:
	enum class File_type
	{
		TEXT,
		BINARY
	};

public:
	explicit File_reader(const std::string& filename, File_type type);

protected:
	template<typename T>
	void read(T&) const;

	template<typename T>
	void read(T*, std::size_t n_elements) const;

protected:
	mutable std::ifstream file_;
};

//////////////////////////////////////////////////////////////////////////

inline File_reader::File_reader(const std::string& filename, File_type type)
{
	if (type == File_type::BINARY)
		file_.open(filename, std::ifstream::binary);
	else
		file_.open(filename);

	if (!file_)
		throw std::runtime_error("File '" + filename + "' not found");
}

template<typename T>
void File_reader::read(T& x) const
{
	file_.read(reinterpret_cast<char*>(&x), sizeof(T));
}

template<typename T>
void File_reader::read(T* x, std::size_t n_elements) const
{
	file_.read(reinterpret_cast<char*>(x), sizeof(T) * n_elements);
}
