//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_map_h
#define load_map_h

#include <iomanip>
#include <iostream>
#include <string>

#include "Decimal.h"

#include "GenMap.hpp"
#include "LoadMap.hpp"


// Genetic map from file

inline Gen::Map load_map(const std::string & filename)
{
	std::cout << "Loading genetic map" << std::endl;
	std::clog << "Loading genetic map" << std::endl;
	
	std::cout << "<< " << filename << std::endl;
	std::clog << "<< " << filename << std::endl;
	
	size_t warn = 0; // count warnings
	
	Gen::Map map;
	
	try
	{
		LoadMap load(filename);
		
		while (load.next())
		{
			try
			{
				Gen::Map::chromosome_t chr;
				Gen::Map::position_t   pos;
				double rate;
				double dist;
				
				if (load.parse(chr, pos, rate, dist))
				{
					map.set(chr, pos, rate, dist);
				}
			}
			catch (const std::string & warning)
			{
				++warn;
				std::clog << "Warning (" << warn << "): " << warning << std::endl;
			}
		}
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	if (warn > 0)
	{
		std::cout << "[" << warn << " warnings - please check log file]" << std::endl;
	}
	
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	return map;
}


// Constant recombination rate

inline Gen::Map load_map(const decimal_t & rr)
{
	std::cout << "No genetic map provided" << std::endl;
	std::clog << "No genetic map provided" << std::endl;
	
	std::cout << "Calculating genetic distances from constant recombination rate: " << std::setprecision(6) << rr << std::endl;
	std::clog << "Calculating genetic distances from constant recombination rate: " << std::setprecision(6) << rr << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	return Gen::Map(rr);
}


#endif /* load_map_h */

