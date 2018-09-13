//
//  Reader.hpp
//  ship
//
//  Created by pkalbers on 14/07/2016.
//  Copyright © 2016 Patrick K. Albers. All rights reserved.
//

#ifndef Reader_hpp
#define Reader_hpp

// C headers (sorted)
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

// C++ headers (sorted)
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>


//
// Read file line by line
// File can be gzip compressed, ending in *.gz
//
class Reader
{
public:
	
	// Line or Token class
	class Current
	{
	private:
		
		friend class Reader;
		
		// construct
		Current(char *, const size_t);
		
		char  *ptr, *beg, *end; // line pointer in buffer, split pointers
		size_t num; // split counter
		
		static constexpr auto split_sep = " \t"; // split at whitespace
		
	public:
		
		// cast
		operator char * () const;
		
		// convert to string
		std::string str() const;
		
		// return line/field size
		size_t size() const;
		
		// return field count
		size_t count() const;
		
		// split line into fields
		bool split(const std::string & = split_sep);
		
		// return field pointer class
		Current field() const;
		
		// convert field to type
		template< typename CAST_T >
		CAST_T convert() const;
		
		const size_t number; // line number
	};
	
	
	// constructs
	Reader(); // initiate, but no filename known yet
	Reader(const std::string &); // open file directly
	Reader(const Reader &) = delete; // do not copy
	Reader & operator = (const Reader &) = delete; // do not assign
	
	// destruct
	~Reader(); // calls close()
	
	// open stream
	void open(const std::string &);
	
	// forward to next line
	bool next();
	
	// reset stream
	void reset();
	
	// close stream
	void close();
	
	// return line pointer class
	Current line() const;
	
	// return line count
	size_t count() const;
	
private:
	
	union Stream
	{
		FILE * file_ptr; // file stream of uncompressed, i.e. text file
		gzFile file_zip; // file stream of compressed, i.e. binary file
	};
	
	typedef std::vector< char >      buffer_t;
	typedef std::vector< char* >     breaks_t;
	typedef breaks_t::const_iterator breaks_p;
	
	// detect gzip compression
	void detect_source();
	
	// read line into buffer
	bool buffer_line();
	
	std::string source; // name of file
	Stream      stream; // file/gzip or string stream
	
	buffer_t    buffer; // stream buffer
	const int   buffer_len; // buffer read length
	char      * buffer_ptr,
			  * buffer_end;
	
	breaks_t    breaks; // vector of newline+1 positions
	breaks_p    breaks_ptr,
	            breaks_end;
	
	size_t num; // line count
	
	bool is_eof, // flag that file is at its end
	     is_open, // flag that stream is open
	     is_compressed; // flag that file is gzip compressed
};



// convert field to type specialisations

template< typename CAST_T >
CAST_T Reader::Current::convert() const
{
	throw std::invalid_argument("Failed conversion in field " + std::to_string(this->num));
	return CAST_T();
}


template<>
inline std::string Reader::Current::convert< std::string >() const
{
	return this->str();
}


template<>
inline int Reader::Current::convert< int >() const
{
	char * end;
	long value = strtol(this->ptr, &end, 10);
	
	if (*end != '\0')
	{
		throw std::invalid_argument("Failed conversion <int> in field " + std::to_string(this->num) + ", token: " + std::string(this->ptr));
	}
	
	return int(value);
}


template<>
inline long Reader::Current::convert< long >() const
{
	char * end;
	long value = strtol(this->ptr, &end, 10);
	
	if (*end != '\0')
	{
		throw std::invalid_argument("Failed conversion <long> in field " + std::to_string(this->num) + ", token: " + std::string(this->ptr));
	}
	
	return value;
}


template<>
inline size_t Reader::Current::convert< size_t >() const
{
	char * end;
	size_t value = strtoul(this->ptr, &end, 10);
	
	if (*end != '\0')
	{
		throw std::invalid_argument("Failed conversion <size_t> in field " + std::to_string(this->num) + ", token: " + std::string(this->ptr));
	}
	
	return value;
}


template<>
inline float Reader::Current::convert< float >() const
{
	char * end;
	float value = strtof(this->ptr, &end);
	
	if (*end != '\0')
	{
		throw std::invalid_argument("Failed conversion <float> in field " + std::to_string(this->num) + ", token: " + std::string(this->ptr));
	}
	
	return value;
}


template<>
inline double Reader::Current::convert< double >() const
{
	char * end;
	double value = strtod(this->ptr, &end);
	
	if (*end != '\0')
	{
		throw std::invalid_argument("Failed conversion <double> in field " + std::to_string(this->num) + ", token: " + std::string(this->ptr));
	}
	
	return value;
}

template<>
inline long double Reader::Current::convert< long double >() const
{
	char * end;
	long double value = strtold(this->ptr, &end);
	
	if (*end != '\0')
	{
		throw std::invalid_argument("Failed conversion <long double> in field " + std::to_string(this->num));
	}
	
	return value;
}


#endif /* Reader_hpp */
