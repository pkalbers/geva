//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_hmm_h
#define load_hmm_h

#include <fstream>
#include <iostream>
#include <string>

#include "Clock.hpp"

#include "IBD.hpp"
#include "IBD_HMM.hpp"
#include "LoadHMM.hpp"


inline IBD::HMM::Model::Data load_hmm(const Gen::Grid::Data grid, const std::string & inits_file, const std::string & emiss_file, const size_t & Ne, const std::string & outfile, const bool print = false)
{
	std::cout << "Generating HMM probabilities from input files" << std::endl;
	std::clog << "Generating HMM probabilities from input files"  << std::endl;
	
	std::cout << "<< " << inits_file << std::endl;
	std::clog << "<< " << inits_file << std::endl;
	
	std::cout << "<< " << emiss_file << std::endl;
	std::clog << "<< " << emiss_file << std::endl;
	
	std::cout << " Effective population size,  Ne: " << Ne << std::endl;
	std::clog << " Effective population size,  Ne: " << Ne << std::endl;
	
	
	const size_t Ns = grid->sample_size() * 2; // number of haplotypes
	
	
	try
	{
		LoadHMM load(Ne, Ns);
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		// initial
		
		std::cout << " Initial state probabilities ... " << std::flush;
		std::clog << " Initial state probabilities ... " << std::flush;
		
		load.make_initial(grid, inits_file);
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		
		// emission
		
		std::cout << " Emission probabilities ... " << std::flush;
		std::clog << " Emission probabilities ... " << std::flush;
		
		load.make_emission(grid, emiss_file);
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		
		// transition
		
		std::cout << " Transition probabilities ... " << std::flush;
		std::clog << " Transition probabilities ... " << std::flush;
		
		load.make_transition(grid);
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		if (print)
		{
			std::cout << "Writing HMM probabilities to file" << std::endl;
			std::clog << "Writing HMM probabilities to file" << std::endl;
			
			const std::string out_inits = outfile + ".initial.txt";
			const std::string out_emiss = outfile + ".emission.txt";
			const std::string out_trans = outfile + ".transition.txt";
		
			std::cout << ">> " << out_inits << std::endl;
			std::clog << ">> " << out_inits << std::endl;
			load.print_initial(out_inits);
			
			std::cout << ">> " << out_emiss << std::endl;
			std::clog << ">> " << out_emiss << std::endl;
			load.print_emission(out_emiss);
			
			std::cout << ">> " << out_trans << std::endl;
			std::clog << ">> " << out_trans << std::endl;
			load.print_transition(out_trans);
			
			std::cout << std::endl;
			std::clog << std::endl;
		}
		
				
		return load.make_model();
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	return nullptr;
}


inline IBD::HMM::Model::Data load_hmm(const Gen::Grid::Data grid, const size_t & Ne)
{
	std::cout << "Generating expected HMM probabilities" << std::endl;
	std::clog << "Generating expected HMM probabilities" << std::endl;
	
	std::cout << " Effective population size,  Ne: " << Ne << std::endl;
	std::clog << " Effective population size,  Ne: " << Ne << std::endl;
	
	const size_t Ns = grid->sample_size() * 2; // number of haplotypes
	
	try
	{
		LoadHMM load(Ne, Ns);
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		// initial
		
		std::cout << " Initial state probabilities ... " << std::flush;
		std::clog << " Initial state probabilities ... " << std::flush;
		
		load.make_initial(grid);
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		
		// emission
		
		std::cout << " Emission probabilities ... " << std::flush;
		std::clog << " Emission probabilities ... " << std::flush;
		
		load.make_emission(grid);
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		
		// transition
		
		std::cout << " Transition probabilities ... " << std::flush;
		std::clog << " Transition probabilities ... " << std::flush;
		
		load.make_transition(grid);
		
		std::cout << "OK" << std::endl;
		std::clog << "OK" << std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		return load.make_model();
	}
	catch (const std::exception & error)
	{
		std::cout << std::endl << "Error: " << error.what() << std::endl;
		std::cerr << error.what() << std::endl << std::endl;
		
		throw std::runtime_error("[Terminated]");
	}
	
	return nullptr;
}


#endif /* load_hmm_h */

