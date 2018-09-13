//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef join_bin_h
#define join_bin_h

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "Gen.hpp"
#include "GenGrid.hpp"


inline std::string join_bin(const std::string & outfile, const std::vector< std::string > & filenames)
{
	std::cout << "Loading data files" << std::endl;
	std::clog << "Loading data files" << std::endl;
	
	const std::string grid_file = outfile + ".bin"; // filename
	
	
	try
	{
		std::vector< std::string >::const_iterator file, file_end = filenames.cend();
		
		// check number of input files
		if (filenames.size() < 2)
			throw std::invalid_argument("No input files to join");
		
		// check non-match between input files and generated output file
		for (file = filenames.cbegin(); file != file_end; ++file)
		{
			if (*file == grid_file)
				throw std::invalid_argument("Output file would have duplicate filename: " + grid_file);
		}
		
		
		Gen::Grid::Join grid(grid_file);
		
		
		// walkabout input files
		for (file = filenames.cbegin(); file != file_end; ++file)
		{
			std::cout << "<< " << *file << " ... " << std::flush;
			std::clog << "<< " << *file << " ... " << std::flush;
			
			Gen::Grid::Data join = std::make_shared<Gen::Grid>(*file);
			
			if (join->marker_size() == 0)
			{
				std::cout << "[empty, skipped]" << std::endl;
				std::clog << "[empty, skipped]" << std::endl;
				continue;
			}
			
			std::cout << "# variants: " << join->marker_size() <<
			//((join->compressed()) ? " (compressed)": "") <<
			std::endl;
			std::clog << "# variants: " << join->marker_size() <<
			//((join->compressed()) ? " (compressed)": "") <<
			std::endl;
			
			grid.insert(join);
		}
		
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		std::cout << "Joining data ... " << std::flush;
		std::clog << "Joining data ... " << std::flush;
		
		grid.finish();
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	return grid_file;
}


#endif /* join_bin_h */

