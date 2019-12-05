#pragma once
#include "file_reader.hpp"
#include <cstddef>
#include <string>
#include <stdexcept>

class Outcar_reader : public File_reader
{
public:
	explicit Outcar_reader(const std::string& filename);

	const std::string& vasp_info() const;
	double fermi_energy() const;

private:
 	void read_file();

private:
	std::string vasp_info_;
	double fermi_energy_;
};

class Bad_outcar_file : public std::runtime_error
{
public:
	Bad_outcar_file(const std::string& err)
		: std::runtime_error("Bad OUTCAR file: " + err)
	{ }
};
