//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "LoadMap.hpp"


// construct

LoadMap::LoadMap(const std::string & filename)
: fields(0)
, count(0)
{
	this->stream.open(filename);

	// read first line
	if (this->stream.next())
	{
		Reader::Current line = this->stream.line();

		// check if header is present (no numeric chars present)
		for (size_t i = 0, end = line.size(); i < end; ++i)
		{
			if (isdigit(line[i]))
			{
				throw std::invalid_argument("Invalid header in genetic map file: " + filename);
			}
		}

		// count number of fields
		while (line.split())
		{
			++this->fields;
		}

		if (this->fields < 3)
		{
			throw std::invalid_argument("Unknown genetic map format: " + filename);
		}
	}
	else
	{
		throw std::invalid_argument("Unable to read from genetic map file: " + filename);
	}
}


// return marker size

size_t LoadMap::size() const
{
	return this->count;
}


// forward to next line

bool LoadMap::next()
{
	return this->stream.next();
}


// parse current line

bool LoadMap::parse(int & chr, size_t & pos, double & rate, double & dist)
{
	if (this->fields == 3) return this->parse_3_fields(chr, pos, rate, dist);
	if (this->fields >= 4) return this->parse_4_fields(chr, pos, rate, dist);
	return false;
}


// parse line with 3 or 4 fields

bool LoadMap::parse_3_fields(int & chr, size_t & pos, double & rate, double & dist)
{
	Reader::Current line = this->stream.line();

	chr  = -1;
	pos  =  0;
	rate = -1;
	dist = -1;

	size_t last_pos = 0;

	while (line.split())
	{
		Reader::Current field = line.field();

		switch (field.number)
		{
			case 0: // position
			{
				pos = field.convert< int >();

				if (pos <= last_pos)
				{
					throw std::invalid_argument("Invalid position on line " + std::to_string(line.number));
				}
				last_pos = pos;

				break;
			}
			case 1: // recombination rate
			{
				rate = field.convert< double >();

				if (rate < 0)
				{
					throw std::invalid_argument("Invalid recombination rate on line " + std::to_string(line.number));
				}

				break;
			}
			case 2: // genetic distance
			{
				dist = field.convert< double >();

				if (dist < 0)
				{
					throw std::invalid_argument("Invalid genetic distance on line " + std::to_string(line.number));
				}

				break;
			}
		}
	}

	return true;
}

bool LoadMap::parse_4_fields(int & chr, size_t & pos, double & rate, double & dist)
{
	Reader::Current line = this->stream.line();

	chr  = -1;
	pos  =  0;
	rate = -1;
	dist = -1;

	size_t last_pos = 0;

	while (line.split())
	{
		Reader::Current field = line.field();

		switch (field.number)
		{
			case 0: // chromosome
			{
				if (field.size() > 3 && field[0] == 'c' && field[1] == 'h' && field[2] == 'r')
				{
					if (field[3] == 'X')
					{
						chr = 23;
					}
					else
					{
						chr = atoi(&field[3]);
					}
				}
				else
				{
					bool flag = false;
					for (size_t i = 0, end = field.size(); i < end; ++i)
					{
						if (isalpha(field[i]))
						{
							flag = true;
							break;
						}
					}

					if (flag)
					{
						throw std::invalid_argument("Invalid chromosome on line " + std::to_string(line.number));
					}

					chr = field.convert< int >();
				}

				break;
			}
			case 1: // position
			{
				pos = field.convert< size_t >();

				if (pos <= last_pos)
				{
					throw std::invalid_argument("Invalid position on line " + std::to_string(line.number));
				}
				last_pos = pos;

				break;
			}
			case 2: // recombination rate
			{
				rate = field.convert< double >();

				if (rate < 0)
				{
					throw std::invalid_argument("Invalid recombination rate on line " + std::to_string(line.number));
				}

				break;
			}
			case 3: // genetic distance
			{
				dist = field.convert< double >();

				if (dist < 0)
				{
					throw std::invalid_argument("Invalid genetic distance on line " + std::to_string(line.number));
				}

				break;
			}
			default:
			{
				return true;
				break;
			}
		}
	}

	return true;
}
