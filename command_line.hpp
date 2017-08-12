#pragma once
#include <string>
#include <vector>

class Command_line
{
public:
	Command_line(int argc, char* argv[]);

	bool is_empty() const;

	bool option_exists(const std::string&) const;
	const std::string& get_option(const std::string&) const;
	const std::string& get_option(const std::string&, const std::string& default_value) const;

private:
	std::vector<std::string> tokens_;
};
