//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef detect_share_h
#define detect_share_h

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "Clock.hpp"

#include "GenGrid.hpp"
#include "GenShare.hpp"


inline bool detect_share(Gen::Share::Data share, Gen::Grid::Data grid, const std::vector<size_t> & list)
{
	std::cout << "Generating sharing index" << std::endl;
	std::clog << "Generating sharing index" << std::endl;
	
	
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
	
	
	Gen::Share::target_t k_list(list.begin(), list.end());
	
	if (k_list.empty())
	{
		throw std::invalid_argument("Invalid target list in sharing detection");
	}
	
	
	try
	{
		std::cout << "...\r" << std::flush;
		
		if (share->detect(k_list, grid))
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


inline bool detect_share(Gen::Share::Data share, Gen::Grid::Data grid, const size_t & min, const size_t & max)
{
	const size_t k_min = std::min(min, max);
	const size_t k_max = std::max(min, max);
	
	const size_t k_beg = (k_min < Gen::Share::minimum) ? Gen::Share::minimum: k_min; // minimum fk
	const size_t k_end = (k_max > grid->sample_size()) ? grid->sample_size(): k_max; // maximum fk
	
	if (k_beg >= k_end)
	{
		throw std::invalid_argument("Invalid target list in sharing detection");
	}
	
	std::vector<size_t> list;
	list.reserve((k_end - k_beg) + 1);
	
	for (size_t i = k_beg; i <= k_end; ++i)
	{
		list.push_back(i);
	}
	
	return detect_share(share, grid, list);
}


#endif /* detect_share_h */

