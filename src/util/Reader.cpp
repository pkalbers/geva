//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Reader.hpp"


//#define READ_BUFFER_SIZE 33554432 // 32 Mb
#define READ_BUFFER_SIZE 268435456 // 256 Mb


//
// Line or Token class
//

// construct

Reader::Current::Current(char * _ptr, const size_t _num)
: ptr(_ptr)
, beg(_ptr)
, end(NULL)
, num(0)
, number(_num)
{}


// cast

Reader::Current::operator char * () const
{
	if (this->ptr != this->beg)
	{
		throw std::runtime_error("Cannot return string pointer after splitting string");
	}
	return this->ptr;
}


// convert to string

std::string Reader::Current::str() const
{
	if (this->ptr != this->beg)
	{
		throw std::runtime_error("Cannot return string after splitting string");
	}
	return std::string(this->ptr);
}


// return line/field size

size_t Reader::Current::size() const
{
	return strlen(this->ptr);
}


// return field count

size_t Reader::Current::count() const
{
	return this->num;
}


// split line into fields

bool Reader::Current::split(const std::string & sep)
{
	// skip leading delimiters
	if (this->beg != NULL)
	{
		while (*this->beg != '\0' && strchr(sep.c_str(), *this->beg) != NULL)
		{
			++this->beg;
		}
		
		if (*this->beg == '\0')
		{
			return false;
		}
	}
	
	if (this->beg == NULL)
	{
		return false;
	}
	
	this->ptr = this->beg; // point to begin
	
	this->end = strpbrk(this->beg, sep.c_str()); // detect delimiter
	
	if (this->end != NULL)
	{
		*this->end = '\0'; // terminate field
		
		this->beg = this->end + 1; // set begin
	}
	else
	{
		this->beg = NULL; // indicate exit
	}
	
	++this->num;
	
	return true;
}


// return field pointer

Reader::Current Reader::Current::field() const
{
	if (this->num == 0)
	{
		throw std::runtime_error("Cannot return field before splitting string");
	}
	return Current(this->ptr, this->num - 1);
}



//
// Read file line by line
// File can be gzip compressed, ending in *.gz
//

// constructs

Reader::Reader()
: buffer(std::min(INT_MAX, READ_BUFFER_SIZE))
, buffer_len(std::min(INT_MAX, READ_BUFFER_SIZE) - 1)
, buffer_ptr(&buffer[0])
, buffer_end(NULL)
, breaks(1)
, breaks_ptr(breaks.cbegin())
, breaks_end(breaks.cend())
, num(0)
, is_eof(false)
, is_open(false)
, is_compressed(false)
{}

Reader::Reader(const std::string & filename)
: Reader()
{
	this->open(filename);
}


// destruct

Reader::~Reader()
{
	this->close();
}


// open stream

void Reader::open(const std::string & filename)
{
	if (this->is_open)
	{
		throw std::runtime_error("File stream already open");
	}
	
	this->source = filename;
	
	// determine file type (text or compressed/binary)
	this->detect_source();
	
	// open file
	if (this->is_compressed)
	{
		if ((this->stream.file_zip = gzopen(this->source.c_str(), "rb")) == NULL)
		{
			throw std::runtime_error("Cannot open compressed file: " + this->source);
		}
	}
	else
	{
		if ((this->stream.file_ptr = fopen(this->source.c_str(), "r")) == NULL)
		{
			throw std::runtime_error("Cannot open file: " + this->source);
		}
	}
	
	this->is_open = true;
}


// forward to next line

bool Reader::next()
{
	if (! this->is_open)
	{
		throw std::runtime_error("Read stream not open");
	}
	
	++this->breaks_ptr;
	
	if (this->breaks_ptr == this->breaks_end)
	{
		this->breaks.clear();
		
		while (this->buffer_line())
		{
			if (this->breaks.size() > 0)
			{
				this->breaks_ptr = this->breaks.cbegin();
				this->breaks_end = this->breaks.cend();
				
				break;
			}
		}
		
		if (this->breaks.size() > 0)
		{
			this->breaks_ptr = this->breaks.cbegin();
			this->breaks_end = this->breaks.cend();
		}
		else
		{
			return false;
		}
	}
	
	++this->num;
	
	return true;
}


// reset stream

void Reader::reset()
{
	if (!this->is_open)
	{
		return;
	}
	
	if (this->is_compressed)
	{
		if (gzrewind(this->stream.file_zip) != 0)
		{
			throw std::runtime_error("Exception while handling compressed file: " + this->source);
		}
	}
	else
	{
		if (fseek(this->stream.file_ptr, 0, SEEK_SET) != 0)
		{
			throw std::runtime_error("Exception while handling file: " + this->source);
		}
	}
	
	this->num = 0;
	this->breaks = std::vector<char*>(1);
	this->breaks_ptr = this->breaks.cbegin();
	this->breaks_end = this->breaks.cend();
}


// close stream

void Reader::close()
{
	if (this->is_open)
	{
		if (this->is_compressed)
			gzclose(this->stream.file_zip);
		else
			fclose(this->stream.file_ptr);
		
		this->is_open = false;
	}
	
	this->buffer.clear();
	this->breaks.clear();
	this->breaks_ptr = this->breaks.cend();
	this->breaks_end = this->breaks.cend();
}


// return line pointer

Reader::Current Reader::line() const
{
	if (this->num == 0)
	{
		throw std::runtime_error("Cannot return line before reading file");
	}
	return Current(*this->breaks_ptr, this->num - 1);
}


// return line count

size_t Reader::count() const
{
	return this->num;
}


// detect gzip compression

void Reader::detect_source()
{
	int c0, c1;
	
	FILE * file = fopen(this->source.c_str(), "rb");
	if (file == NULL)
	{
		throw std::runtime_error("Cannot open file: " + this->source);
	}
	
	c0 = fgetc(file);
	c1 = fgetc(file);
	
	if (feof(file) || ferror(file))
	{
		throw std::runtime_error("Cannot read file: " + this->source);
	}
	
	// magic number
	this->is_compressed = (c0 == 0x1f && c1 == 0x8b);
	
	fclose(file);
}


// read line into buffer

bool Reader::buffer_line()
{
	if (this->is_eof)
	{
		return false;
	}
	
	int n = 0;
	const size_t size = strlen(this->buffer_ptr);
	
	// copy overhang to front of buffer
	if (size > 0 && this->buffer_ptr != &this->buffer[0])
	{
		memmove(&this->buffer[0], this->buffer_ptr, size);
	}
	
	// increase buffer size to fit overhang and next read
	if (size + this->buffer_len > this->buffer.size() - 1)
	{
		this->buffer.resize(size + this->buffer_len + 1);
	}
	
	this->buffer_ptr = &this->buffer[size];
	*this->buffer_ptr = '\0';
	
	if (this->is_compressed)
	{
		// read from gzip
		n = gzread(this->stream.file_zip,
				   this->buffer_ptr,
				   this->buffer_len);
	}
	else
	{
		// read from file
		n = (int)fread(this->buffer_ptr,
					   sizeof(this->buffer[0]),
					   this->buffer_len,
					   this->stream.file_ptr);
	}
	
	// handle end of file
	if (n < 1)
	{
		if (strlen(&this->buffer[0]) > 0)
		{
			this->breaks.push_back(&this->buffer[0]);
		}
		
		this->is_eof = true;
		
		return false;
	}
	
	*(this->buffer_ptr + n) = '\0'; // terminate read
	
	this->buffer_ptr = &this->buffer[0];
	this->buffer_end = strchr(this->buffer_ptr, '\n'); // detect newline
	
	// walkabout buffer
	while (this->buffer_end != NULL)
	{
		this->breaks.push_back(this->buffer_ptr);
		
		*this->buffer_end = '\0'; // replace newline
		
		this->buffer_ptr = this->buffer_end + 1;
		this->buffer_end = strchr(this->buffer_ptr, '\n'); // detect newline
	}
	
	return true;
}

