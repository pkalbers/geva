//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef infer_age_h
#define infer_age_h

#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Progress.hpp"
#include "Clock.hpp"

#include "GenGrid.hpp"
#include "GenShare.hpp"

#include "AgeInfer.hpp"

#include "Threadpool.h"


inline void infer_age(const Age::Param::Data param,
					  const IBD::DetectMethod method,
					  const decimal_t max_miss,
					  const std::string & output,
					  const bool print_cle_full,
					  const bool print_ccf_full,
					  const Gen::Share::Data share,
					  const Gen::Grid::Data grid,
					  Clock & time,
					  const IBD::HMM::Model::Data hmm_model = nullptr,
					  const IBD::SIM::Result::Data simres = nullptr,
					  const size_t threads = 1,
					  const size_t batch_limit = 1000)
{
	std::cout << "Age estimation, using ";
	std::clog << "Age estimation, using ";

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
//			if (hmm_model->do_iterative)
//			{
//				std::cout << "Iterative approach" << std::endl;
//				std::clog << "Iterative approach" << std::endl;
//			}
			break;
		case IBD::DETECT_SIM:
			std::cout << IBD::SIM::name << std::endl;
			std::clog << IBD::SIM::name << std::endl;
			break;
		case IBD::DETECT_VOID: return;
	}


	// print results to files

	const std::string file_pairs = output + ".pairs.txt";
	const std::string file_sites = output + ".sites.txt";
	
	
	// detection method
	/*
	switch (method)
	{
		case IBD::DETECT_DGT:
			file_pairs += '.' + IBD::DGT::abbr;
			file_sites += '.' + IBD::DGT::abbr;
			break;
		case IBD::DETECT_FGT:
			file_pairs += '.' + IBD::FGT::abbr;
			file_sites += '.' + IBD::FGT::abbr;
			break;
		case IBD::DETECT_HMM:
			file_pairs += '.' + IBD::HMM::abbr;
			file_sites += '.' + IBD::HMM::abbr;
			if (hmm_model->do_iterative)
			{
				file_pairs += 'i';
				file_sites += 'i';
			}
			break;
		case IBD::DETECT_SIM:
			file_pairs += '.' + IBD::SIM::abbr;
			file_sites += '.' + IBD::SIM::abbr;
			break;
		case IBD::DETECT_VOID: return;
	}
	*/
	
	// age estimation clocks
	
	//if (!param->use_mut_clock && !param->use_rec_clock)
	//	throw std::invalid_argument("At least one clock must be included to perform age estimation");
	
	//std::string clock_str = "";
	
	//if (param->use_mut_clock) clock_str += "Mut";
	//if (param->use_rec_clock) clock_str += "Rec";
	
	//file_pairs += '.' + clock_str;
	//file_sites += '.' + clock_str;
	/*
	if (param->use_hard_brks)
	{
		file_pairs += ".HardBreaks";
		file_sites += ".HardBreaks";
	}
	else
	{
		file_pairs += ".SoftBreaks";
		file_sites += ".SoftBreaks";
	}
	*/
	
	const std::string file_pairs_distr = output + ".pairs.distr.txt";
	const std::string file_sites_distr = output + ".sites.distr.txt";
	
	//file_pairs += ".txt";
	//file_sites += ".txt";

	std::cout << ">> " << file_pairs << ((print_ccf_full) ? " + " + file_pairs_distr: "") << std::endl;
	std::cout << ">> " << file_sites << ((print_cle_full) ? " + " + file_sites_distr: "") << std::endl;
	
	std::clog << ">> " << file_pairs << ((print_ccf_full) ? " + " + file_pairs_distr: "") << std::endl;
	std::clog << ">> " << file_sites << ((print_cle_full) ? " + " + file_sites_distr: "") << std::endl;

	std::cout << std::endl;
	std::clog << std::endl;


	// print parameters

	std::cout << " Age estimation parameters" << std::endl;
	
	std::cout << "  Effective population size, Ne:  " << static_cast<size_t>(param->Ne) << std::endl;
	std::cout << "  Mutation rate per site:         " << std::scientific << std::setprecision(4) << param->Mr << std::endl;
	std::cout << "  # max. concordant pairs:        " << param->limit_sharers << std::endl;
	std::cout << "  # max. discordant pairs:        " << param->outgroup_size << std::endl;
	
	//std::cout << "   [" << ((param->use_hard_brks) ? "X": " ") << "] Hard breakpoints (assumes correct IBD)" <<  std::endl;
	//std::cout << "   [" << ((param->use_mut_clock) ? "X": " ") << "] Mutation clock" <<  std::endl;
	//std::cout << "   [" << ((param->use_rec_clock) ? "X": " ") << "] Recombination clock" <<  std::endl;
//	std::cout << "   [" << ((param->apply_nearest_neighb) ? "X": " ") << "] Prioritise nearest neighbours" <<  std::endl;
	//std::cout << "   [" << ((param->use_tree_consistency) ? "X": " ") << "] Count only consistent differences (mutation clock)" <<  std::endl;
	//std::cout << "   [" << ((param->apply_filter_fixed) ? "X": " ") << "] Filter pairs: exclude fixed proportion" <<  std::endl;
	//std::cout << "   [" << ((param->apply_filter_detect) ? "X": " ") << "] Filter pairs: automatic threshold detection" <<  std::endl;
	
//	if (method == IBD::DETECT_HMM)
//	{
//		std::cout << "   [" << ((param->use_post_prob) ? "X": " ") << "] Include HMM posterior probability" <<  std::endl;
//	}
	
	//std::cout << "  Data required: " << ((param->use_hard_brks || param->use_mut_clock || method == IBD::DetectMethod::DETECT_FGT) ? "Haplotypes (phased)": "Genotypes (phased or unphased)") << std::endl;

	std::cout << std::endl;
	
	
	//std::cout << " Max. pairwise missing rate: " << std::fixed << std::setprecision(4) << max_miss << std::endl;
	//std::clog << " Max. pairwise missing rate: " << std::fixed << std::setprecision(4) << max_miss << std::endl;
	
	if (threads > 1)
	{
		std::cout << " # threads = " << threads << std::endl;
		std::clog << " # threads = " << threads << std::endl;
	}


	try
	{
		// queue for detection and inference

		Age::Queue queue(share, batch_limit, grid, param, simres);
		
		std::cout << " # pairwise analyses = " << queue.size() << std::endl;
		std::clog << " # pairwise analyses = " << queue.size() << std::endl;
		
		std::cout << std::endl;
		std::clog << std::endl;
		
		
		// write results to streams
		
		std::ofstream stream_pairs(file_pairs);
		std::ofstream stream_sites(file_sites);
		
		Age::Pair::print_header(stream_pairs, param, false);
		Age::Site::print_header(stream_sites, param, false);
		
		
		std::ofstream stream_pairs_distr;
		std::ofstream stream_sites_distr;
		
		if (print_cle_full)
		{
			stream_sites_distr.open(file_sites_distr);
			Age::Site::print_header(stream_sites_distr, param, true);
		}
		
		if (print_ccf_full)
		{
			stream_pairs_distr.open(file_pairs_distr);
			Age::Pair::print_header(stream_pairs_distr, param, true);
		}
		
		
		// walkabout queue in batches
		
		//size_t n_pairs_done = 0;
		//size_t n_pairs_good = 0;
		
		//size_t n_sites_done = 0;
		//size_t n_sites_good = 0;
		
		Progress prog(queue.size());
		
		size_t warn = 0;
		
		while (true)
		{
			// get batch from queue
			
			const size_t n_queue = queue.next(simres);
			
			if (n_queue == 0)
				break;
			
			
			// detect
			
			Age::Pair::List::iterator pair, pair_end = queue.pairs.end();
			
			if (threads > 1)
			{
				Threadpool< Age::Infer > pool(threads, &Age::Infer::run);
				
				size_t count = 0;
				
				for (pair = queue.pairs.begin(); pair != pair_end; ++pair)
				{
					pool.task(Age::Infer(param, method, max_miss, *pair, grid, hmm_model));
					
					if (++count == batch_limit)
						break;
				}
				
				pool.open(&prog, &time); // start other threads
				pool.exec(&prog); // execute on this thread
				pool.wait(); // wait for completion of all threads
			}
			else
			{
				size_t count = 0;
				
				size_t x_con = 0;
				size_t x_dis = 0;
				
				for (pair = queue.pairs.begin(); pair != pair_end; ++pair)
				{
					if ((*pair)->sharing)
						++x_con;
					else
						++x_dis;
					
					Age::Infer infer(param, method, max_miss, *pair, grid, hmm_model);
					
					prog.update();
					
					infer.run();
					
					if (++count == batch_limit)
						break;
				}				
			}
			
			
			// infer
			
			Age::Site::List::iterator site, site_end = queue.sites.end();
			
			for (site = queue.sites.begin(); site != site_end; ++site)
			{
				try
				{
					(*site)->estimate(param);
				}
				catch (const std::string & warning)
				{
					std::cerr << "Warning: " << warning << std::endl;
					++warn;
				}
			}
			
			
			// print sites
			
			for (site = queue.sites.begin(); site != site_end; ++site)
			{
				if ((*site)->done)
				{
					(*site)->print(stream_sites, param->Ne, false, false); // raw
					(*site)->print(stream_sites, param->Ne, false, true);  // adj
					
					if (print_cle_full)
					{
						(*site)->print(stream_sites_distr, param->Ne, true, false); // raw
						(*site)->print(stream_sites_distr, param->Ne, true, true);  // adj
					}
					
					
					// print pairs
					
					Age::Pair::List::iterator p, p_end = (*site)->list.end();
					
					for (p = (*site)->list.begin(); p != p_end; ++p)
					{
						(*p)->print(stream_pairs, param, false);
						
						if (print_ccf_full)
							(*p)->print(stream_pairs_distr, param, true);
					}
				}
			}
		}
		
		
		prog.finish();
		
		//const double f_pairs = static_cast<double>(100) * static_cast<double>(n_pairs_good) / static_cast<double>(n_pairs_done);
		//const double f_sites = static_cast<double>(100) * static_cast<double>(n_sites_good) / static_cast<double>(n_sites_done);
		
		//std::cout << "# sucessful pairwise analyses = " << n_pairs_good << " of " << n_pairs_done << " (" << std::setprecision(1) << f_pairs << "%)" << std::endl;
		//std::clog << "# sucessful pairwise analyses = " << n_pairs_good << " of " << n_pairs_done << " (" << std::setprecision(1) << f_pairs << "%)" << std::endl;
		
		//std::cout << "# completed focal sites = " << n_sites_good << " of " << n_sites_done << " (" << std::setprecision(1) << f_sites << "%)" << std::endl;
		//std::clog << "# completed focal sites = " << n_sites_good << " of " << n_sites_done << " (" << std::setprecision(1) << f_sites << "%)" << std::endl;
		
		if (warn > 0)
		{
			std::cout << "[" << warn << " warnings - please check error file]" << std::endl;
		}
		
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


#endif /* infer_age_h */

