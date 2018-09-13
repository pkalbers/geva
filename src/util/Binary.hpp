//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Binary_hpp
#define Binary_hpp

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>


// Read and write binary files
class Binary
{
private:
	
	static_assert(CHAR_BIT * sizeof(float)  == 32, "System not compatible (float != 32 bits)");
	static_assert(CHAR_BIT * sizeof(double) == 64, "System not compatible (double != 64 bits)");
	
	// specialisation for allowed (portable) types
	template<typename T> struct Constrain : std::false_type {};
	
	
	typedef std::unique_ptr<std::FILE, int (*)(std::FILE *)> file_ptr;
	typedef std::unordered_map< size_t, std::fpos_t >        file_map;
	
	// Make file with custom deleter for writing + updating
	static file_ptr make(const std::string & filename);
	
	// Open file with custom deleter for reading + updating
	static file_ptr open(const std::string & filename);
	
	const std::string name; // filename
	file_ptr          file; // file pointer
	file_map           map; // indexed file positions
	
	const bool auto_delete;
	
public:
	
	enum mode { WRITE, READ };
	
	// constructs
	Binary(const std::string &, const mode, const bool = false);
	Binary();
	Binary(Binary &&); // move
	Binary(const Binary &) = delete; // no copy
	
	// destruct
	~Binary();
	
	// get filename
	std::string filename() const;
	
	
	// Write
	
	// Write value to file
	template< typename WriteAs, typename IsType, typename std::enable_if<std::is_arithmetic<IsType>::value, int>::type = 0 >
	void write(const IsType value)
	{
		static_assert(Constrain<WriteAs>::value, "Write type not allowed");
		
		if (std::is_same<WriteAs, IsType>::value)
		{
			if (std::fwrite(&value, sizeof(WriteAs), 1, this->file.get()) == 1)
				return;
		}
		else
		{
			const WriteAs _value = static_cast<WriteAs>(value);
			if (std::fwrite(&_value, sizeof(WriteAs), 1, this->file.get()) == 1)
				return;
		}
		
		throw std::runtime_error("Error while writing value to binary file");
	}
	
	// Write array to file
	template< typename WriteAs, typename IsType, typename std::enable_if<std::is_arithmetic<IsType>::value, int>::type = 0 >
	void write(const IsType * data_ptr, const size_t length = 1)
	{
		static_assert(Constrain<WriteAs>::value, "Write type not allowed");
		
		if (std::is_same<WriteAs, IsType>::value)
		{
			if (std::fwrite(data_ptr, sizeof(WriteAs), length, this->file.get()) == length)
				return;
		}
		else
		{
			WriteAs buffer[length];
			for (size_t i = 0; i < length; ++i)
			{
				buffer[i] = static_cast<WriteAs>(data_ptr[i]);
			}
			if (std::fwrite(buffer, sizeof(WriteAs), length, this->file.get()) == length)
				return;
		}
		
		throw std::runtime_error("Error while writing array to binary file");
	}
	
	
	// Read
	
	// Read value from file, direct assign
	template< typename ReadAs, typename ToType, typename std::enable_if<std::is_arithmetic<ToType>::value, int>::type = 0 >
	ToType read()
	{
		static_assert(Constrain<ReadAs>::value, "Read type not allowed");
		
		ReadAs value;
		if (std::fread(&value, sizeof(ReadAs), 1, this->file.get()) != 1)
		{
			throw std::runtime_error("Error while reading value from binary file");
		}
		return (std::is_same<ReadAs, ToType>::value) ? value: static_cast<ToType>(value);
	}
	
	// Read value from file
	template< typename ReadAs, typename ToType, typename std::enable_if<std::is_arithmetic<ToType>::value, int>::type = 0 >
	void read(ToType & value)
	{
		static_assert(Constrain<ReadAs>::value, "Read type not allowed");
		
		if (std::is_same<ReadAs, ToType>::value)
		{
			if (std::fread(&value, sizeof(ReadAs), 1, this->file.get()) != 1)
			{
				throw std::runtime_error("Error while reading value from binary file");
			}
		}
		else
		{
			value = read<ReadAs, ToType>();
		}
	}
	
	// Read array from file
	template< typename ReadAs, typename ToType, typename std::enable_if<std::is_arithmetic<ToType>::value, int>::type = 0 >
	void read(ToType * data_ptr, const size_t length = 1)
	{
		static_assert(Constrain<ReadAs>::value, "Read type not allowed");
		
		if (std::is_same<ReadAs, ToType>::value)
		{
			if (std::fread(data_ptr, sizeof(ReadAs), length, this->file.get()) != length)
			{
				throw std::runtime_error("Error while reading array from binary file");
			}
		}
		else
		{
			ReadAs buffer[length];
			if (std::fread(buffer, sizeof(ReadAs), length, this->file.get()) != length)
			{
				throw std::runtime_error("Error while reading array from binary file");
			}
			for (size_t i = 0; i < length; ++i)
			{
				data_ptr[i] = static_cast<ToType>(buffer[i]);
			}
		}
	}
	
	
	// Match
	
	// Read and match to value
	template< typename ReadAs, typename ToType, typename std::enable_if<std::is_arithmetic<ToType>::value, int>::type = 0 >
	void match(const ToType value)
	{
		const ToType found = read<ReadAs, ToType>();
		if (found != value)
		{
			throw std::runtime_error("False alignment of value in binary file");
		}
	}
	
	// Read and match to array
	template< typename ReadAs, typename ToType, typename std::enable_if<std::is_arithmetic<ToType>::value, int>::type = 0 >
	void match(const ToType * data_ptr, const size_t length = 1)
	{
		ToType buffer[length];
		read<ReadAs, ToType>(buffer, length);
		for (size_t i = 0; i < length; ++i)
		{
			if (buffer[i] != data_ptr[i])
			{
				throw std::runtime_error("False alignment of array in binary file");
			}
		}
	}
	
	
	// Navigate
	
	// Seek file position
	template< typename SkipType >
	void skip(const size_t length = 1)
	{
		static_assert(Constrain<SkipType>::value, "Read type not allowed");
		
		if (sizeof(SkipType) * static_cast<unsigned long long>(length) > LONG_MAX)
		{
			throw std::runtime_error("Integral type cannot address file position");
		}
		
		if (std::fseek(this->file.get(), sizeof(SkipType) * length, SEEK_CUR) != 0)
		{
			throw std::runtime_error("Error while navigating binary file position");
		}
	}
	
	// Get current file position
	void here(const file_map::key_type);
	
	// Go to file position
	void jump(const file_map::key_type);
	
	// Roll back to begin of file
	void reset();
};


// specialisation for allowed (portable) types

template<> struct Binary::Constrain<bool> : std::true_type {};
template<> struct Binary::Constrain<char> : std::true_type {};

template<> struct Binary::Constrain<int8_t>  : std::true_type {};
template<> struct Binary::Constrain<int16_t> : std::true_type {};
template<> struct Binary::Constrain<int32_t> : std::true_type {};
template<> struct Binary::Constrain<int64_t> : std::true_type {};

template<> struct Binary::Constrain<uint8_t>  : std::true_type {};
template<> struct Binary::Constrain<uint16_t> : std::true_type {};
template<> struct Binary::Constrain<uint32_t> : std::true_type {};
template<> struct Binary::Constrain<uint64_t> : std::true_type {};

template<> struct Binary::Constrain<float>  : std::true_type {};
template<> struct Binary::Constrain<double> : std::true_type {};


#endif /* Binary_hpp */

