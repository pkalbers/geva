//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Command_hpp
#define Command_hpp

#include <array>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
#include <utility>




namespace Command
{
	// Parameter types
	template < typename > class Value;
	template < typename > class Vector;
	template < typename, size_t > class Array;
	template < typename, typename... > class Tuple;
	class Bool;
	
	
	static constexpr size_t print_width = 10;
	
	
	// Line parser
	class Line
	{
	public:
		
		// construct
		Line(int, const char **);
		
		// parse command line
		void parse();
		
		// print command line
		void print(std::ostream & = std::cout) const;
		
		// clean up parsing
		void finish() const;
		
		
		// get parameter types
		
		template < typename T >
		bool get(Value<T> &, const bool = false, const T = T());
		
		template < typename T >
		bool get(Vector<T> &, const bool = false, const std::vector<T> = std::vector<T>());
		
		template < typename T, size_t N >
		bool get(Array<T, N> &, const bool = false, const std::array<T, N> = std::array<T, N>());
		
		template < typename T0, typename T1 >
		bool get(Tuple<T0, T1> &, const bool = false, const std::tuple<T0, T1> = std::tuple<T0, T1>());
		
		template < typename T0, typename T1, typename T2 >
		bool get(Tuple<T0, T1, T2> &, const bool = false, const std::tuple<T0, T1, T2> = std::tuple<T0, T1, T2>());
		
		template < typename T0, typename T1, typename T2, typename T3 >
		bool get(Tuple<T0, T1, T2, T3> &, const bool = false, const std::tuple<T0, T1, T2, T3> = std::tuple<T0, T1, T2, T3>());
		
		template < typename T0, typename T1, typename T2, typename T3, typename T4 >
		bool get(Tuple<T0, T1, T2, T3, T4> &, const bool = false, const std::tuple<T0, T1, T2, T3, T4> = std::tuple<T0, T1, T2, T3, T4>());
		
		bool get(Bool &, const bool = false, const bool = bool());
		
		
	private:
		
		static constexpr char flag = '-';
		
		using line_t = std::vector< std::string >;
		using pars_t = std::vector< line_t::const_pointer >;
		using args_t = std::unordered_map< char, pars_t >;
		using opts_t = std::unordered_map< std::string, pars_t >;
		
		enum Mode : int
		{
			PAR,
			ARG,
			OPT
		};
		
		
		// check if argument was provided
		bool match(const char, const std::string &) const;
		
		// retrieve argument parameter
		pars_t fetch(const char, const std::string &, const int = -1);
		
		
		// convert to type
		
		template <typename T>
		static T convert(line_t::const_pointer ptr, std::true_type)
		{
			double value;
			std::istringstream iss(*ptr);
			iss >> value;
			if (iss.fail())
				throw std::runtime_error(std::string("Unable to interpret command line parameter: ") + iss.str());
			return static_cast<T>(value);
		}
		
		template <typename T>
		static T convert(line_t::const_pointer ptr, std::false_type)
		{
			T value;
			std::istringstream iss(*ptr);
			iss >> value;
			if (iss.fail())
				throw std::runtime_error(std::string("Unable to interpret command line parameter: ") + iss.str());
			return value;
		}
		
		template <typename T>
		static T convert(line_t::const_pointer ptr)
		{
			return convert<T>(ptr, std::is_arithmetic<T>());
		}
		
		
		// make string of parameter defintion
		static std::string defstr(const char, const std::string &);
		
		
		const std::string path; // execution path
		const line_t      line; // parsed command line
		
		pars_t pars;
		args_t args;
		opts_t opts;
		
		
		template < typename > friend class Value;
		template < typename > friend class Vector;
		template < typename, size_t > friend class Array;
		template < typename, typename... > friend class Tuple;
		friend class Bool;
		
		static bool good;
		static bool help;
	};
	
	
	//
	// Parameter types
	//
	
	// Value
	
	template < typename T >
	class Value
	{
		friend class Line;
		
		template < typename > friend class Vector;
		template < typename, size_t > friend class Array;
		template < typename, typename... > friend class Tuple;
		friend class Bool;
		
		
	private:
		
		typedef T type;
		
		// print help
		void if_help() const
		{
			if (Line::help && this->valid)
			{
				std::cout << std::left << std::setw(print_width + 5) << Line::defstr(this->arg, this->opt) << " : ";
				std::cout << this->txt << std::endl << std::endl;
			}
		}
		
		const char        arg;
		const std::string opt;
		const std::string txt;
		
		bool valid;
		bool input;
		
		
	public:
		
		T value;
		
		
		// constructs
		Value(char && argument, std::string && option, std::string && description)
		: arg(argument)
		, opt(option)
		, txt(description)
		, valid(false)
		, input(false)
		{
			// check if definition is valid
			if (isalpha(this->arg) && !this->opt.empty())
			{
				for (size_t i = 0; i < this->opt.size(); ++i)
				{
					if (isalnum(this->opt[i]))
						break;
				}
				this->valid = true;
			}
			this->if_help();
		}
		Value(char && argument, std::string && description)
		: arg(argument)
		, opt(std::string())
		, txt(description)
		, valid(false)
		, input(false)
		{
			if (isalpha(this->arg))
			{
				this->valid = true;
			}
			this->if_help();
		}
		Value(std::string && option, std::string && description)
		: arg('\0')
		, opt(option)
		, txt(description)
		, valid(false)
		, input(false)
		{
			if (!this->opt.empty())
			{
				for (size_t i = 0; i < this->opt.size(); ++i)
				{
					if (isalnum(this->opt[i]))
						break;
				}
				this->valid = true;
			}
			this->if_help();
		}
		Value(std::string && description)
		: arg('\0')
		, opt(std::string())
		, txt(description)
		, valid(true)
		, input(false)
		{
			this->if_help();
		}
		
		// casts
		operator T const & () const
		{
			return this->value;
		}
		T const & operator () () const
		{
			return this->value;
		}
		
		// assigns
		T & operator = (const T & other)
		{
			this->value = other;
			return this->value;
		}
		T & operator = (T && other)
		{
			this->value = std::move(other);
			return this->value;
		}
		
		// check if provided
		bool good() const { return this->input; }
	};
	
	
	// Vector
	
	template < typename T >
	class Vector : public Value< std::vector< T > >
	{
		using self = Value< std::vector< T > >;
		
	public:
		
		// construct
		Vector(char && argument, std::string && option, std::string && description)
		: self(std::move(argument), std::move(option), std::move(description))
		{}
		Vector(char && argument, std::string && description)
		: self(std::move(argument), std::move(description))
		{}
		Vector(std::string && option, std::string && description)
		: self(std::move(option), std::move(description))
		{}
		Vector(std::string && description)
		: self(std::move(description))
		{}
		
		// cast
		T const & operator [] (const size_t i) const
		{
			return this->value.at(i);
		}
		
		// return size
		size_t size() const
		{
			return this->value.size();
		}
	};
	
	
	// Array
	
	template < typename T, size_t N >
	class Array : public Value< std::array< T, N > >
	{
		using self = Value< std::array< T, N > >;
		
	public:
		
		// construct
		Array(char && argument, std::string && option, std::string && description)
		: self(std::move(argument), std::move(option), std::move(description))
		{}
		Array(char && argument, std::string && description)
		: self(std::move(argument), std::move(description))
		{}
		Array(std::string && option, std::string && description)
		: self(std::move(option), std::move(description))
		{}
		Array(std::string && description)
		: self(std::move(description))
		{}
		
		// cast
		constexpr T const & operator [] (const size_t i) const
		{
			return this->value[i];
		}
	};
	
	
	// Tuple
	
	template < typename T, typename ...Ts >
	class Tuple : public Value< std::tuple< T,  Ts... > >
	{
		using self = Value< std::tuple< T,  Ts... > >;
		
	public:
		
		// construct
		Tuple(char && argument, std::string && option, std::string && description)
		: self(std::move(argument), std::move(option), std::move(description))
		{}
		Tuple(char && argument, std::string && description)
		: self(std::move(argument), std::move(description))
		{}
		Tuple(std::string && option, std::string && description)
		: self(std::move(option), std::move(description))
		{}
		Tuple(std::string && description)
		: self(std::move(description))
		{}
		
		// cast
		template < size_t N >
		constexpr typename std::tuple_element< N, std::tuple< Ts... > >::type const & get() const
		{
			return std::get<N>(this->value);
		}
	};
	
	
	// Bool
	
	class Bool : public Value< bool >
	{
		using self = Value< bool >;
		
	public:
		
		// construct
		Bool(char && argument, std::string && option, std::string && description)
		: self(std::move(argument), std::move(option), std::move(description))
		{
			this->value = false;
		}
		Bool(char && argument, std::string && description)
		: self(std::move(argument), std::move(description))
		{
			this->value = false;
		}
		Bool(std::string && option, std::string && description)
		: self(std::move(option), std::move(description))
		{
			this->value = false;
		}
		Bool(std::string && description)
		: self(std::move(description))
		{
			this->value = false;
		}
	};
	
	
	
	//
	// get parameter types
	//
	
	template < typename T >
	bool Line::get(Value<T> & param, const bool req, const T def)
	{
		static_assert(std::is_convertible<int,float>::value, "Default value is not convertible to defined type");
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, 1);
			
			param.value = Line::convert<T>(input.front());
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	
	template < typename T >
	bool Line::get(Vector<T> & param, const bool req, const std::vector<T> def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt);
			
			pars_t::const_iterator it, ti = input.cend();
			for (it = input.cbegin(); it != ti; ++it)
			{
				param.value.push_back(Line::convert<T>(*it));
			}
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	
	template < typename T, size_t N >
	bool Line::get(Array<T, N> & param, const bool req, const std::array<T, N> def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, N);
			
			for (size_t i = 0, end = input.size(); i != end; ++i)
			{
				param.value[i] = Line::convert<T>(input[i]);
			}
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	template < typename T0, typename T1 >
	bool Line::get(Tuple<T0, T1> & param, const bool req, const std::tuple<T0, T1> def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, 2);
			
			std::get<0>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1> > >(input[0]);
			std::get<1>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1> > >(input[1]);
			
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	template < typename T0, typename T1, typename T2 >
	bool Line::get(Tuple<T0, T1, T2> & param, const bool req, const std::tuple<T0, T1, T2> def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, 2);
			
			std::get<0>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2> > >(input[0]);
			std::get<1>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2> > >(input[1]);
			std::get<2>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2> > >(input[2]);
			
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	template < typename T0, typename T1, typename T2, typename T3 >
	bool Line::get(Tuple<T0, T1, T2, T3> & param, const bool req, const std::tuple<T0, T1, T2, T3> def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, 2);
			
			std::get<0>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3> > >(input[0]);
			std::get<1>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3> > >(input[1]);
			std::get<2>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3> > >(input[2]);
			std::get<3>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3> > >(input[3]);
			
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	template < typename T0, typename T1, typename T2, typename T3, typename T4 >
	bool Line::get(Tuple<T0, T1, T2, T3, T4> & param, const bool req, const std::tuple<T0, T1, T2, T3, T4> def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, 2);
			
			std::get<0>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3, T4> > >(input[0]);
			std::get<1>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3, T4> > >(input[1]);
			std::get<2>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3, T4> > >(input[2]);
			std::get<3>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3, T4> > >(input[3]);
			std::get<4>(param.value) = Line::convert< std::tuple_element< 0, std::tuple<T0, T1, T2, T3, T4> > >(input[4]);
			
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	
	
	inline
	bool Line::get(Bool & param, const bool req, const bool def)
	{
		if (!Line::good)  throw std::logic_error("Command line was not parsed");
		if (!param.valid) throw std::invalid_argument("Invalid format in parameter definition:" + Line::defstr(param.arg, param.opt));
		if (param.input)  throw std::logic_error("Command line parameter already fetched:" + Line::defstr(param.arg, param.opt));
		
		if (this->match(param.arg, param.opt))
		{
			pars_t input = this->fetch(param.arg, param.opt, 0);
			
			if (input.size() > 0)
			{
				throw std::runtime_error("Unexpected input after command line flag: " + Line::defstr(param.arg, param.opt));
			}
			
			param.value = true;
			param.input = true;
			return true;
		}
		
		if (req) throw std::runtime_error("Command line parameter required:" + Line::defstr(param.arg, param.opt));
		
		param.value = def;
		param.input = false;
		return false;
	}
	

	
}


#endif /* Command_hpp */

