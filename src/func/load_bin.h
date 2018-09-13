//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_bin_h
#define load_bin_h


#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "Gen.hpp"
#include "GenGrid.hpp"


inline Gen::Grid::Data load_bin(const std::string & filename)
{
	std::cout << "Loading data" << std::endl << "<< " << filename << std::endl;
	std::clog << "Loading data" << std::endl << "<< " << filename << std::endl;
	
	try
	{
		Gen::Grid::Data grid = std::make_shared<Gen::Grid>(filename);
		
		
		size_t n_phased = 0;
		
		Gen::Sample::Iterator it, ti = grid->sample().cend();
		
		for (it = grid->sample().cbegin(); it != ti; ++it)
		{
			if (it->phase)
			{
				++n_phased;
			}
		}
		
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		std::cout << " # individuals: " << grid->sample().size() << " (" << n_phased << " phased)" << std::endl;
		std::cout << " # variants:    " << grid->marker().size() <<
		//((grid->compressed()) ? " (compressed)": "") <<
		std::endl;
		
		std::clog << " # individuals: " << grid->sample().size() << " (" << n_phased << " phased)" << std::endl;
		std::clog << " # variants:    " << grid->marker().size() <<
		//((grid->compressed()) ? " (compressed)": "") <<
		std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		return grid;
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	return nullptr;
}


#endif /* load_bin_h */
