//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef ibd_by_site_h
#define ibd_by_site_h

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "Progress.hpp"
#include "Clock.hpp"

#include "GenGrid.hpp"
#include "GenShare.hpp"

#include "IBDbySite.hpp"

#include "Threadpool.h"


inline IBD::Site::Result::Data ibd_by_site(const IBD::DetectMethod method, const decimal_t max_miss, const std::string & output, const Gen::Share::Data share, const Gen::Grid::Data grid, Clock & time, const IBD::HMM::Model::Data hmm_model = nullptr, const size_t threads = 1)
{
	std::cout << "IBD estimation, using ";
	std::clog << "IBD estimation, using ";
	
	switch (method)
	{
		case IBD::DETECT_DGT:
			std::cout << IBD::DGT::name << std::endl;
			std::clog << IBD::DGT::name << std::endl;
			break;
		case IBD::DETECT_FGT:
			std::cout << IBD::FGT::name << std::endl;
			std::clog << IBD::FGT::name << std::endl;
			break;
		case IBD::DETECT_HMM:
			std::cout << IBD::HMM::name << std::endl;
			std::clog << IBD::HMM::name << std::endl;
			if (hmm_model->do_iterative)
			{
				std::cout << "Iterative approach" << std::endl;
				std::clog << "Iterative approach" << std::endl;
			}
			break;
		case IBD::DETECT_SIM:
			throw std::logic_error("But why?");
			break;
		case IBD::DETECT_VOID: return nullptr;
	}
	
	std::cout << " Max. pairwise missing rate = " << std::setprecision(4) << max_miss << std::endl;
	std::clog << " Max. pairwise missing rate = " << std::setprecision(4) << max_miss << std::endl;
	
	if (threads > 1)
	{
		std::cout << " # threads = " << threads << std::endl;
		std::clog << " # threads = " << threads << std::endl;
	}
	
	
	try
	{
		// prepare result container
		
		IBD::Site::Result::Data result = std::make_shared< IBD::Site::Result >(share);
		
		std::cout << " # scans = " << result->size() << std::endl;
		std::clog << " # scans = " << result->size() << std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		// execute
		
		Progress prog(result->get().size(), "sites");
		
		IBD::Site::Target::List::const_iterator list, end = result->get().cend();
		
		if (threads > 1)
		{
			Threadpool< IBD::Site::Detect > pool(threads - 1, &IBD::Site::Detect::run);
			
			for (list = result->get().cbegin(); list != end; ++list)
			{
				pool.task(IBD::Site::Detect(method, max_miss, *list, grid, hmm_model));
			}
			
			pool.open(&prog, &time); // start other threads
			pool.exec(&prog); // execute on this thread
			pool.wait(); // wait for completion of all threads
		}
		else
		{
			for (list = result->get().cbegin(); list != end; ++list)
			{
				IBD::Site::Detect detect(method, max_miss, *list, grid, hmm_model);
				
				prog.update();
				
				detect.run();
			}
		}
		
		prog.finish();
		time.print(std::cout);
		time.print(std::clog);
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		// print results to file
		
		std::string filename = output + ".ibd.";
		
		switch (method)
		{
			case IBD::DETECT_DGT: filename += IBD::DGT::abbr + ".txt"; break;
			case IBD::DETECT_FGT: filename += IBD::FGT::abbr + ".txt"; break;
			case IBD::DETECT_HMM: if (hmm_model->do_iterative) { filename += IBD::HMM::abbr + "i.txt"; } else { filename += IBD::HMM::abbr + ".txt"; } break;
			case IBD::DETECT_SIM: throw std::logic_error("This cannot happen!"); break;
			case IBD::DETECT_VOID: return nullptr;
		}
		
		std::cout << "Writing results to file" << std::endl << ">> " << filename << std::endl;
		std::clog << "Writing results to file" << std::endl << ">> " << filename << std::endl;
		
		std::ofstream file(filename);
		
		const size_t lines = result->print(file);
		
		std::cout << "# lines written = " << lines << std::endl;
		std::clog << "# lines written = " << lines << std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		return result;
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	
	return nullptr;
}


#endif /* ibd_by_site_h */

