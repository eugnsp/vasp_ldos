#include "command_line.hpp"
#include <cstddef>
#include <algorithm>
#include <stdexcept>

Command_line::Command_line(int argc, char* argv[])
{
	for (int i = 1; i < argc; ++i)
		tokens_.push_back(argv[i]);
}

bool Command_line::is_empty() const
{
	return tokens_.empty();
}

bool Command_line::option_exists(const std::string& option) const
{
	return std::find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
}

const std::string& Command_line::get_option(const std::string& option) const
{
	auto pos = std::find(tokens_.begin(), tokens_.end(), option);
	if (pos != tokens_.end() && ++pos != tokens_.end())
		return *pos;

	throw std::runtime_error("Option '" + option + "' does not exist");
}

const std::string& Command_line::get_option(
	const std::string& option,
	const std::string& default_value) const
{
	auto pos = std::find(tokens_.begin(), tokens_.end(), option);
	if (pos != tokens_.end() && ++pos != tokens_.end())
		return *pos;
	else
		return default_value;
}