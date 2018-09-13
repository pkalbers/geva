//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef LoadMap_hpp
#define LoadMap_hpp

#include <ctype.h>
#include <stdio.h>

#include <exception>
#include <utility>
#include <set>
#include <sstream>

#include <iostream>

#include "Reader.hpp"


class LoadMap
{
public:
	
	// construct
	LoadMap(const std::string &);
	
	// return marker size
	size_t size() const;
	
	// forward to next line
	bool next();
	
	// parse current line
	bool parse(int &, size_t &, double &, double &);
	
private:
	
	// parse line with 3 or 4 fields
	bool parse_3_fields(int &, size_t &, double &, double &);
	bool parse_4_fields(int &, size_t &, double &, double &);
	
	Reader stream; // stream of Map file
	size_t fields; // number of detected fields
	size_t count;  // marker counter
};


#endif /* LoadMap_hpp */
