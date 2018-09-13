//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_share_h
#define load_share_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "Clock.hpp"

#include "GenGrid.hpp"
#include "GenShare.hpp"


// save sharing index

inline void save_share(const Gen::Share::Data share, const std::string & prefix, const Gen::Grid::Data grid)
{
	std::cout << "Saving sharing index" << std::endl;
	std::clog << "Saving sharing index" << std::endl;
	
	const std::string filename = prefix + ".sharing.idx";
	
	std::cout << ">> " << filename << std::endl;
	std::clog << ">> " << filename << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	try
	{
		share->save(filename, grid);
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
}


// load sharing index

inline bool load_share(const std::string & filename, Gen::Share::Data share, const Gen::Grid::Data grid)
{
	std::cout << "Loading sharing index" << std::endl;
	std::clog << "Loading sharing index" << std::endl;
	
	std::cout << "<< " << filename << std::endl;
	std::clog << "<< " << filename << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	if (share->region.first != 0 || share->region.second != 0)
	{
		const size_t lower = std::min(share->region.first, share->region.second);
		const size_t upper = std::max(share->region.first, share->region.second);
		
		std::cout << " Chromosome region:  Position " << lower << " to " << upper << " (inclusive)" << std::endl;
		std::clog << " Chromosome region:  Position " << lower << " to " << upper << " (inclusive)" << std::endl;
		
		share->region.first  = lower;
		share->region.second = upper;
	}
	
	if (share->max_sites > 0)
	{
		std::cout << " # max. sites: " << share->max_sites << std::endl;
		std::clog << " # max. sites: " << share->max_sites << std::endl;
	}
	
	if (share->max_pairs > 0)
	{
		std::cout << " # max. pairs: " << share->max_pairs << std::endl;
		std::clog << " # max. pairs: " << share->max_pairs << std::endl;
	}
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	try
	{
		if (share->load(filename, grid))
		{
			share->print(grid->sample_size() * 2, std::cout);
			share->print(grid->sample_size() * 2, std::clog);
			
			std::cout << std::endl;
			std::clog << std::endl;
			
			return true;
		}
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	return false;
}


#endif /* load_share_h */
