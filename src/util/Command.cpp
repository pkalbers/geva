//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Command.hpp"


#define COMMAND_HELP_KEYWORD "help"


using namespace Command;


bool Line::good = false;
bool Line::help = false;


// construct

Line::Line(int argc, const char ** argv)
: path(argv[0])
, line(argv + 1, argv + argc)
{
	Line::good = true;
	
	for (size_t i = 0; i < this->line.size(); ++i)
	{
		if (this->line[i] == std::string(1, flag) + COMMAND_HELP_KEYWORD ||
			this->line[i] == std::string(2, flag) + COMMAND_HELP_KEYWORD)
		{
			Line::help = true;
			break;
		}
	}
}


// parse command line

void Line::parse()
{
	// exit on help
	if (Line::help)
	{
		throw int(0);
	}
	
	char        last_arg;
	std::string last_opt;
	
	Mode current = PAR; // leading parameter
	
	for (size_t i = 0; i < this->line.size(); ++i)
	{
		const std::string value = this->line[i];
		
		// argument
		if (value.size() == 2 && value[0] == flag && value[1] != flag && isalpha(value[1]))
		{
			last_arg = value[1];
			
			if (this->args.count(last_arg) != 0)
			{
				throw std::runtime_error(std::string("Duplicate command line argument: ") + value);
			}
			this->args[last_arg] = pars_t();
			
			current = ARG;
			continue;
		}
		
		// option
		if (value.size() > 2 && value[0] == flag && value[1] == flag && value[2] != flag && isalpha(value[2]))
		{
			last_opt = std::string(value.cbegin() + 2, value.cend());
			
			if (this->opts.count(last_opt) != 0)
			{
				throw std::runtime_error(std::string("Duplicate command line option: ") + value);
			}
			this->opts[last_opt] = pars_t();
			
			current = OPT;
			continue;
		}
		
		// allocate line index to value
		switch (current)
		{
			case PAR: this->pars.push_back(&this->line[i]); break;
			case ARG: this->args[last_arg].push_back(&this->line[i]); break;
			case OPT: this->opts[last_opt].push_back(&this->line[i]); break;
		}
	}
}


// print command line

void Line::print(std::ostream & stream) const
{
	stream << "Path: " << this->path << std::endl;
	stream << "Call:";
	for (size_t i = 0; i < this->line.size(); ++i)
	{
		stream << ' ' << this->line[i];
	}
	stream << std::endl << std::endl;
}


// clean up parsing

void Line::finish() const
{
	if (!Line::good) throw std::logic_error("Command line was not parsed");
	
	if (this->pars.size() == 0 && this->args.size() == 0 && this->opts.size() == 0)
	{
		return;
	}
	
	std::ostringstream oss;
	
	// parameters
	if (this->pars.size() > 0)
	{
		oss << std::endl;
		for (pars_t::const_iterator it = this->pars.cbegin(), end = this->pars.cend(); it != end; ++it)
		{
			oss << ' ' << **it;
		}
	}
	
	// arguments
	if (this->args.size() > 0)
	{
		oss << std::endl;
		for (args_t::const_iterator it = this->args.cbegin(), end = this->args.cend(); it != end; ++it)
		{
			oss << ' ' << flag << it->first;
		}
	}
	
	// options
	if (this->opts.size() > 0)
	{
		oss << std::endl;
		for (opts_t::const_iterator it = this->opts.cbegin(), end = this->opts.cend(); it != end; ++it)
		{
			oss << ' ' << flag << flag << it->first;
		}
	}
	
	throw std::runtime_error(std::string("Unused command line arguments: ") + oss.str());
}


// check if argument was provided

bool Line::match(const char arg, const std::string & opt) const
{
	bool is_arg = false;
	bool is_opt = false;
	
	if (arg == '\0' && opt.empty() && this->pars.size() > 0)
	{
		return true;
	}
	
	if (arg != '\0' && this->args.count(arg) != 0)
	{
		is_arg = true;
	}
	
	if (!opt.empty() && this->opts.count(opt) != 0)
	{
		is_opt = true;
	}
	
	if (is_arg && is_opt)
	{
		throw std::runtime_error(std::string("Redundant defintion of command line arguments: ") + "-" + arg + " and --" + opt);
	}
	
	if (is_arg || is_opt)
	{
		return true;
	}
	
	return false;
}


// retrieve argument parameter

Line::pars_t Line::fetch(const char arg, const std::string & opt, const int limit)
{
	pars_t input;
	
	if (arg == '\0' && opt.empty())
	{
		this->pars.swap(input);
		this->pars.clear();
	}
	
	if (this->args.count(arg) != 0)
	{
		this->args[arg].swap(input);
		this->args.erase(arg);
	}
	
	if (this->opts.count(opt) != 0)
	{
		this->opts[opt].swap(input);
		this->opts.erase(opt);
	}
	
	if (limit > -1 && static_cast<int>(input.size()) != limit)
	{
		throw std::runtime_error(std::string("Unexpected defintion of command line parameter: ") + defstr(arg, opt) + "\n"
								 "Expected number of parameters: " + std::to_string(limit) + "\n"
								 "Detected number of paratemers: " + std::to_string(input.size()));
	}
	
	return input;
}


// make string of parameter defintion

std::string Line::defstr(const char arg, const std::string & opt)
{
	const std::string arg_str = (arg == '\0') ? std::string(): std::string(" ") + std::string(1, Line::flag) + arg;
	const std::string opt_str = (opt.empty()) ? std::string(): std::string(" ") + std::string(2, Line::flag) + opt;
	
	if (arg_str.empty() && opt_str.empty())
	{
		return " [leading parameter] ";
	}
	
	return arg_str + opt_str;
}


