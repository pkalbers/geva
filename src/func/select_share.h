//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef select_share_h
#define select_share_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "Clock.hpp"

#include "GenGrid.hpp"
#include "GenShare.hpp"


inline bool select_share(Gen::Share::Data share, Gen::Grid::Data grid, const std::string & filename)
{
	std::cout << "Scanning target positions" << std::endl;
	std::clog << "Scanning target positions" << std::endl;
	
	std::cout << "<< " << filename << std::endl;
	std::clog << "<< " << filename << std::endl;
	
	
//	if (share->region.first != 0 || share->region.second != 0)
//	{
//		const size_t lower = std::min(share->region.first, share->region.second);
//		const size_t upper = std::max(share->region.first, share->region.second);
//
//		std::cout << " Chromosome region:  Position " << lower << " to " << upper << " (inclusive)" << std::endl;
//		std::clog << " Chromosome region:  Position " << lower << " to " << upper << " (inclusive)" << std::endl;
//
//		share->region.first  = lower;
//		share->region.second = upper;
//	}
//
//	if (share->max_sites > 0)
//	{
//		std::cout << " # max. sites: " << share->max_sites << std::endl;
//		std::clog << " # max. sites: " << share->max_sites << std::endl;
//	}
//
//	if (share->max_pairs > 0)
//	{
//		std::cout << " # max. pairs: " << share->max_pairs << std::endl;
//		std::clog << " # max. pairs: " << share->max_pairs << std::endl;
//	}
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	try
	{
		// read input file
		
		std::ifstream input(filename, std::ifstream::in);
		
		std::string str;
		size_t      pos;
		
		std::set<size_t> positions; // list of target sites
		
		while (input >> str)
		{
			std::istringstream stream(str);
			
			if (stream >> pos)
			{
				const size_t prev = positions.size();
				
				positions.insert(pos);
				
				if (positions.size() == prev)
				{
					throw std::invalid_argument("Duplicate positions in input file");
				}
				
				continue;
			}
			
			throw std::invalid_argument("Unable to convert '" + str + "' to integral type");
		}
		
		const size_t n_provided = positions.size();
		
		if (n_provided == 0)
		{
			std::runtime_error("No target sites in input file");
		}
		
		
		std::cout << " # target sites: " << n_provided << std::endl;
		std::clog << " # target sites: " << n_provided << std::endl;
		
		
		std::cout << "...\r" << std::flush;
		
		// create share table
		const size_t n_matched = share->select(positions, grid);
		
		
		if (n_provided != n_matched)
		{
			std::cout << " [" << (n_provided - n_matched) << " positions not viable]" << std::endl;
			std::clog << " [" << (n_provided - n_matched) << " positions not viable]" << std::endl;
		}
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		// print summary and return
		
		if (n_matched > 0)
		{
//			share->print(grid->sample_size() * 2, std::cout);
//			share->print(grid->sample_size() * 2, std::clog);
//
//			std::cout << std::endl;
//			std::clog << std::endl;
			
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


inline bool select_share(Gen::Share::Data share, Gen::Grid::Data grid, const size_t & pos)
{
	std::cout << "Scanning target position: " << pos << std::endl;
	std::clog << "Scanning target position: " << pos << std::endl;
	
	
//	if (share->region.first != 0 || share->region.second != 0)
//	{
//		const size_t lower = std::min(share->region.first, share->region.second);
//		const size_t upper = std::max(share->region.first, share->region.second);
//
//		std::cout << " Chromosome region:  Position " << lower << " to " << upper << " (inclusive)" << std::endl;
//		std::clog << " Chromosome region:  Position " << lower << " to " << upper << " (inclusive)" << std::endl;
//
//		share->region.first  = lower;
//		share->region.second = upper;
//	}
//
//	if (share->max_sites > 0)
//	{
//		std::cout << " # max. sites: " << share->max_sites << std::endl;
//		std::clog << " # max. sites: " << share->max_sites << std::endl;
//	}
//
//	if (share->max_pairs > 0)
//	{
//		std::cout << " # max. pairs: " << share->max_pairs << std::endl;
//		std::clog << " # max. pairs: " << share->max_pairs << std::endl;
//	}
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	try
	{
		std::set<size_t> position;
		
		position.insert(pos);
		
		
		std::cout << "...\r" << std::flush;
		
		// create share table
		const size_t matched = share->select(position, grid);
		
		
		if (matched == 0)
		{
			std::cout << " [Position not viable]" << std::endl;
			std::clog << " [Position not viable]" << std::endl;
		}
		else
		{
//			share->print(grid->sample_size() * 2, std::cout);
//			share->print(grid->sample_size() * 2, std::clog);
//
//			std::cout << std::endl;
//			std::clog << std::endl;
			
			return true;
		}
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	return false;
}


#endif /* select_share_h */

