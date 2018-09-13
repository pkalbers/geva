//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Binary.hpp"


// constructs

Binary::Binary(const std::string & filename, const Binary::mode m, const bool remove)
try
: name(filename)
, file((m == Binary::WRITE) ? Binary::make(filename): Binary::open(filename))
, auto_delete(remove)
{}
catch (const std::exception & ex)
{
	throw std::runtime_error(std::string("Error while opening binary file: ") + filename + std::string("\n") + ex.what());
}

Binary::Binary()
try
: name(std::string())
, file(std::tmpfile(), std::fclose)
, auto_delete(false)
{}
catch (const std::exception & ex)
{
	throw std::runtime_error(std::string("Error while creating temporary file\n") + ex.what());
}

Binary::Binary(Binary && other)
try
: name(std::move(other.name))
, file(std::move(other.file))
, map(std::move(other.map))
, auto_delete(other.auto_delete)
{}
catch (const std::exception & ex)
{
	throw std::runtime_error(std::string("Error while handling binary file\n") + ex.what());
}


// destruct

Binary::~Binary()
{
	if (this->auto_delete && this->file)
	{
		this->file.reset();
		std::remove(this->name.data());
	}
}


// get filename

std::string Binary::filename() const
{
	return this->name;
}


// Make file with custom deleter for writing + updating

Binary::file_ptr Binary::make(const std::string & filename)
{
	return file_ptr(std::fopen(filename.data(), "wb+"), std::fclose);
}


// Open file with custom deleter for reading + updating

Binary::file_ptr Binary::open(const std::string & filename)
{
	return file_ptr(std::fopen(filename.data(), "rb"), std::fclose);
}


// Get current file position

void Binary::here(const Binary::file_map::key_type index)
{
	if (this->map.count(index) != 0)
	{
		throw std::runtime_error("Binary map index already in use: " + std::to_string(index));
	}
	
	std::fpos_t pos;
	if (std::fgetpos(this->file.get(), &pos) != 0)
	{
		throw std::runtime_error("Error while retrieving binary file position");
	}
	this->map[index] = pos;
}


// Go to file position

void Binary::jump(const Binary::file_map::key_type index)
{
	if (this->map.count(index) == 0)
	{
		throw std::runtime_error("Unknown binary map index: " + std::to_string(index));
	}
	
	if (std::fsetpos(this->file.get(), &this->map[index]) != 0)
	{
		throw std::runtime_error("Error while locating binary file position");
	}
}


// Roll back to begin of file

void Binary::reset()
{
	rewind(this->file.get());
}



















