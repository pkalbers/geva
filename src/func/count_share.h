//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef count_share_h
#define count_share_h

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Progress.hpp"

#include "GenGrid.hpp"
#include "GenVariant.hpp"


inline void count_share(const std::string & output, const Gen::Grid::Data grid, const std::vector<size_t> & list)
{
	std::cout << "Generating pairwise sharing matrix" << std::endl;
	std::clog << "Generating pairwise sharing matrix" << std::endl;
	
	const std::string prefix = output + ".psm.";
	
	std::cout << ">> " << prefix << "*.txt" << std::endl;
	std::clog << ">> " << prefix << "*.txt" << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	try
	{
		const size_t n_pairs = (grid->sample_size() * (grid->sample_size() - 1)) / 2;
		
		
		// walkabout list
		
		std::vector<size_t>::const_iterator fk, fk_end = list.cend();
		
		for (fk = list.cbegin(); fk != fk_end; ++fk)
		{
			const size_t k = *fk;
			const std::string name = prefix + std::to_string(k) + ".txt";
			
			
			// prepare index vector
			
			std::vector<size_t> index;
			index.reserve(grid->marker_size());
			
			Gen::Marker::Vector::const_iterator marker, marker_end = grid->marker().cend();
			
			for (marker = grid->marker().cbegin(); marker != marker_end; ++marker)
			{
				if (marker->hap_count[Gen::H1] == k)
				{
					index.push_back(marker->index.value);
				}
			}
			
			
			// output file
			
			std::ofstream file(name);
			
			file << "SampleID0 SampleID1 N" << std::endl; // print header
			
			
			std::cout << " # shared: " << *fk << std::endl;
			std::clog << " # shared: " << *fk << " ... " << std::flush;
			
			Progress prog(n_pairs);
			
			
			// count for each pair
			
			for (size_t a = 0; a < grid->sample_size() - 1; ++a)
			{
				const Gen::Variant::Vector::Data va = grid->get(a); // read
				
				for (size_t b = a + 1; b < grid->sample_size(); ++b)
				{
					const Gen::Variant::Vector::Data vb = grid->get(b); // read
					
					size_t count = 0;
					
					prog.update();
					
					
					std::vector<size_t>::const_iterator idx, end = index.cend();
					
					for (idx = index.cbegin(); idx != end; ++idx)
					{
						if (Gen::is_genotype<Gen::G1>(va->gen(*idx)) &&
							Gen::is_genotype<Gen::G1>(vb->gen(*idx)))
							++count;
					}
					
					
					file << a << ' ' << b << ' ' << count << std::endl; // print
				}
			}
			
			
			file.close();
			
			prog.finish();
			
			std::cout << std::endl;
			std::clog << "OK" << std::endl;
		}
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
}


inline void count_share(const std::string & output, const Gen::Grid::Data grid, const size_t & min, const size_t & max)
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
	
	return count_share(output, grid, list);
}


#endif /* count_share_h */

