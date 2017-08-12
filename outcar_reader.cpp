#include "outcar_reader.hpp"
#include <cctype>
#include <string>
#include <limits>
#include <algorithm>

namespace
{
void l_trim(std::string& s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
		[](int ch) { return !std::isspace(ch); }));
}

// trim from end (in place)
void r_trim(std::string& s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(),
		[](int ch) { return !std::isspace(ch); }).base(), s.end());
}

bool begins_with_ci(const std::string& s, const std::string& match)
{
	return s.size() >= match.size() &&
		std::equal(match.begin(), match.end(), s.begin(),
			[](int ch1, int ch2) { return std::tolower(ch1) == std::tolower(ch2); });
}
}

Outcar_reader::Outcar_reader(const std::string& filename)
	: File_reader(filename, File_type::TEXT)
{	
	read_file();
}

const std::string& Outcar_reader::vasp_info() const
{
	return vasp_info_;
}

double Outcar_reader::fermi_energy() const
{
	return fermi_energy_;
}

void Outcar_reader::read_file()
{
	std::getline(file_, vasp_info_);
	l_trim(vasp_info_);
	r_trim(vasp_info_);
	if (!begins_with_ci(vasp_info_, "vasp"))
		throw Bad_outcar_file("Bad header");

	std::string line;
	while (std::getline(file_, line))
	{
		l_trim(line);
		if (begins_with_ci(line, "e-fermi"))
		{
			auto colon = std::find(line.begin(), line.end(), ':');
			if (colon == line.end() || ++colon == line.end())
				continue;

			try
			{
				fermi_energy_ = std::stod(std::string{colon, line.end()});
			}
			catch (...)
			{
				fermi_energy_ = std::numeric_limits<double>::quiet_NaN();
			}
			break;
		}
	}
}
