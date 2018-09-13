//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_gen_h
#define load_gen_h

#include <iostream>
#include <stdexcept>
#include <string>

#include "Clock.hpp"
#include "Progress.hpp"
#include "LoadGen.hpp"
#include "Gen.hpp"
#include "GenGrid.hpp"


inline std::string load_gen(const std::string & gen_file, const std::string & sample_file, const Gen::Map & gmap, const std::string & outfile,
							const size_t & buffer_limit, const bool compress, const bool local_tmp = false, const size_t chunk_beg = 0, const size_t chunk_end = 0)
{
	std::cout << "Loading genotype data from GEN + SAMPLE files to binary file" << std::endl;
	std::clog << "Loading genotype data from GEN + SAMPLE files to binary file" << std::endl;
	
	std::cout << "<< " << gen_file << std::endl;
	std::clog << "<< " << gen_file << std::endl;
	
	std::cout << "<< " << sample_file << std::endl;
	std::clog << "<< " << sample_file << std::endl;
	
	
	const std::string grid_file = outfile + ".bin"; // filename
	
	std::cout << ">> " << grid_file << std::endl;
	std::clog << ">> " << grid_file << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	size_t warn = 0; // count warnings
	
	try
	{
		Clock   time;
		LoadGen load(gen_file, sample_file); // reads sample vector
		
		const size_t sample_size = load.sample.size();
		const size_t buffer_size = static_cast<size_t>(static_cast<double>(sample_size * buffer_limit) / static_cast<double>(1024 * 1024 * sizeof(Gen::value_t)));
		
		std::cout << " Detected sample size: " << sample_size << " individuals" << std::endl;
		std::clog << " Detected sample size: " << sample_size << " individuals" << std::endl;
		
		std::cout << " Buffer window size: " << buffer_limit << " variants (~" << buffer_size << " Mb)" << std::endl;
		std::clog << " Buffer window size: " << buffer_limit << " variants (~" << buffer_size << " Mb)" << std::endl;
		
		std::cout << " Data compression: " << (compress ? "On": "Off") << std::endl;
		std::clog << " Data compression: " << (compress ? "On": "Off") << std::endl;
		
		if (chunk_beg != chunk_end)
		{
			load.filter.position_beg = std::min(chunk_beg, chunk_end);
			load.filter.position_end = std::max(chunk_beg, chunk_end);
			
			std::cout << " Chunk: positions from " << load.filter.position_beg << " (inclusive) to " << load.filter.position_end << " (exclusive)" << std::endl;
			std::clog << " Chunk: positions from " << load.filter.position_beg << " (inclusive) to " << load.filter.position_end << " (exclusive)" << std::endl;
		}
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		// make new grid
		
		Gen::Grid::Make buffer(grid_file, sample_size, buffer_limit, compress);
		
		std::clog << " Running" << std::endl;
		
		Progress progress("variants");
		
		size_t ntmp = 0; // number of temporary files
		
		while (load.next())
		{
			try
			{
				if (load.parse(buffer, gmap))
				{
					progress.update();
					
					if (buffer.full)
					{
						progress.halt();
						
						std::cout << " Saving buffer to temporary file ..." << std::flush;
						std::clog << " Saving buffer to temporary file ..." << std::flush;
						
						buffer.save(local_tmp);
						++ntmp;
						
						std::cout << " OK" << std::endl;
						std::clog << " OK" << std::endl;
					}
				}
			}
			catch (const std::string & warning)
			{
				++warn;
				std::clog << "Warning (" << warn << "): " << warning << std::endl;
			}
		}
		
		
		// concatenate files
		
		progress.halt();
		
		if (ntmp == 0)
		{
			std::cout << " Writing output file ..." << std::flush;
			std::clog << " Writing output file ..." << std::flush;
		}
		else
		{
			std::cout << " Combining temporary files ..." << std::flush;
			std::clog << " Combining temporary files ..." << std::flush;
		}
		
		buffer.finish(load.sample, load.marker);
		
		std::cout << " OK" << std::endl;
		std::clog << " OK" << std::endl;
		
		
		progress.finish();
		std::clog << " Done" << std::endl;
		time.print(std::clog);
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
	
	
	return grid_file;
}


#endif /* load_gen_h */

