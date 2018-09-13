//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef print_ibd_h
#define print_ibd_h

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


inline void print_ibd(const IBD::DetectMethod method, const std::string & output, const Gen::Share::Data share, const Gen::Grid::Data grid, const IBD::HMM::Model::Data hmm_model = nullptr)
{
	std::cout << "Write IBD estimation to file, using ";
	std::clog << "Write IBD estimation to file, using ";
	
	std::string filename = output + ".ibd.distr.";
	
	switch (method)
	{
		case IBD::DETECT_DGT:
			std::cout << IBD::DGT::name << std::endl;
			std::clog << IBD::DGT::name << std::endl;
			filename += IBD::DGT::abbr;
			break;
		case IBD::DETECT_FGT:
			std::cout << IBD::FGT::name << std::endl;
			std::clog << IBD::FGT::name << std::endl;
			filename += IBD::FGT::abbr;
			break;
		case IBD::DETECT_HMM:
			std::cout << IBD::HMM::name << std::endl;
			std::clog << IBD::HMM::name << std::endl;
			filename += IBD::HMM::abbr;
			break;
		case IBD::DETECT_SIM:
			throw std::logic_error("Nope!");
			break;
		case IBD::DETECT_VOID: return;
	}
	
	filename += ".txt";
	
	std::cout << ">> " << filename << std::endl;
	std::clog << ">> " << filename << std::endl;
	
	std::cout << std::endl;
	std::clog << std::endl;
	
	
	std::ofstream file(filename);
	
	// header
	switch (method)
	{
		case IBD::DETECT_DGT: IBD::DGT::Algorithm::print_header(file); break;
		case IBD::DETECT_FGT: IBD::FGT::Algorithm::print_header(file); break;
		case IBD::DETECT_HMM: IBD::HMM::Algorithm::print_header(file); break;
		case IBD::DETECT_SIM: throw std::logic_error("This should never show up!"); break;
		case IBD::DETECT_VOID: return;
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
		
		for (list = result->get().cbegin(); list != end; ++list)
		{
			IBD::Site::Detect detect(method, 0.0, *list, grid, hmm_model);
			
			prog.update();
			
			detect.print(file);
		}
		
		prog.finish();
		
		std::cout << std::endl;
		std::clog << std::endl;
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
}


#endif /* print_ibd_h */

