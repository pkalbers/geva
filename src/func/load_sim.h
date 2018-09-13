//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_sim_h
#define load_sim_h

#include "Progress.hpp"
#include "LoadSim.hpp"


inline IBD::SIM::Result::Data load_sim(const std::string & filename)
{
	std::cout << "Loading IBD results" << std::endl;
	std::clog << "Loading IBD results" << std::endl;
	
	std::cout << "<< " << filename << std::endl;
	std::clog << "<< " << filename << std::endl;
	
	IBD::SIM::Result::Data result = std::make_shared<IBD::SIM::Result>();
	
	Progress progress;
	
	try
	{
		LoadSim load(filename);
		
		while (load.next())
		{
			size_t site;
			IBD::SIM::Truth element;
			
			progress.update();
			
			if (load.parse(site, element))
			{
				result->set(site, element);
				continue;
			}
			
			throw std::runtime_error("Unable to read file");
		}
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	progress.finish();
	
	// count focal sites
	const size_t n = result->get().size();
	
	std::cout << "Number of focal sites: " << n << std::endl;
	std::clog << "Number of focal sites: " << n << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	return result;
}


#endif /* load_sim_h */
