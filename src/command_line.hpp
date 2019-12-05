#pragma once
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

class Command_line
{
public:
	Command_line(int argc, char* argv[])
	{
		for (int i = 1; i < argc; ++i)
			tokens_.push_back(argv[i]);
	}

	bool is_empty() const
	{
		return tokens_.empty();
	}

	bool option_exists(const std::string& option) const
	{
		return std::find(tokens_.begin(), tokens_.end(), option) != tokens_.end();
	}

	const std::string& get_option(const std::string& option) const
	{
		auto pos = std::find(tokens_.begin(), tokens_.end(), option);
		if (pos != tokens_.end() && ++pos != tokens_.end())
			return *pos;

		throw std::runtime_error("Option '" + option + "' does not exist");
	}

	const std::string& get_option_or(const std::string& option, const std::string& default_value) const
	{
		auto pos = std::find(tokens_.begin(), tokens_.end(), option);
		if (pos != tokens_.end() && ++pos != tokens_.end())
			return *pos;
		else
			return default_value;
	}

private:
	std::vector<std::string> tokens_;
};
