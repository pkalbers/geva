//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include <algorithm>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <exception>
#include <stdexcept>

#include "load_map.h"
#include "load_vcf.h"
//#include "load_gen.h"
//#include "load_hap.h"
#include "load_bin.h"
//#include "join_bin.h"
#include "load_hmm.h"
//#include "load_sim.h"

//#include "detect_share.h"
#include "select_share.h"
#include "load_share.h"

//#include "ibd_by_pair.h"
//#include "ibd_by_site.h"
//#include "print_ibd.h"
//#include "count_share.h"

#include "infer_age.h"

#include "Command.hpp"
#include "Redirect.hpp"
#include "Clock.hpp"


//
// Main
//
int main(int argc, const char * argv[])
{
	Command::Line line(argc, argv);
	
	// common arguments
	Command::Value< std::string >  output('o', "out", "Prefix for generated output files");
	Command::Value< size_t >       thread('t', "treads", "Number of threads for parallel execution (default: 1)");
	Command::Value< size_t >       buffer('b', "buffer", "Memory buffer, upper limit in megabytes (default: no limit)");
	Command::Value< size_t >       seed("seed", "Set seed for random operations (may not be reproducible if t>1)");
	
	// input arguments
	Command::Value< std::string >    input_bin_file('i', "input", "Pre-processed binary input file");
	Command::Value< std::string >    input_vcf_file("vcf", "VCF input file (optionally GZIP compressed)");
	Command::Value< std::string >    input_map_file("map", "Gen map input file (optionally GZIP compressed)");
	Command::Value< double >         input_rec_rate("rec", "Recombination rate, per site per generation (default: 1e-08)");
	Command::Value< size_t >         input_lines("maxLines", "Number of lines to be buffered while processing input file (default: 500000)");
	Command::Bool                    local_tmp_files("localTmpFiles", "Temporary files are stored in local directory when parsing input file");
	Command::Array< std::string, 2 > input_hmm_file("hmm", "Hidden Markov Model, 2 input files: (1) empirical initial state and (2) emission probabilties");
	Command::Value< std::string >    input_anc_file("AncAlleles", "File containing chromosome, position, and ancestral allele (no header, tab-separated)");
	Command::Bool                    require_anc_match("removeUnmatchedAncestral", "Remove variants where the ancestral allele does not match ref/alt");
	Command::Value< bool >           require_snp("onlySNP", "Remove non-SNP variants");
	
	// genomic position arguments
	Command::Value< size_t >      share_position("position", "Target position");
	Command::Value< std::string > share_batch("positions", "Batch file containing target positions (any white-space separation)");
	
	// age estimation parameters
	Command::Value< size_t > effective_size("Ne", "Effective population size (default: 10000)");
	Command::Value< double > mutation_rate("mut", "Mutation rate, per site per generation (default: 1e-08)");
	Command::Value< size_t > age_limit_sharers("maxConcordant", "Maximim number of concordant pairs to be selected (default: 100)");
	Command::Value< size_t > age_outgroup_size("maxDiscordant", "Maximum number of discordant pairs to be selected (default: 100)");
	Command::Value< double > max_missing("maxMissing", "Maximum missing proportion of sites in HMM (default: 5%)");
	
	
	// parse command line
	try
	{
		line.parse();
		
		line.get(output, true);
		line.get(thread, false, size_t(1));
		line.get(buffer, false, std::numeric_limits<size_t>::max()); // default: all individuals
		line.get(seed, false);
		
		
		if (!line.get(input_bin_file, false)) // BIN
		{
			// VCF file
			line.get(input_vcf_file, false);
			
			// genetic map
			line.get(input_map_file, false);
			
			// optionally constant recombination rate
			line.get(input_rec_rate, false, double(1e-08)); // default: 1e-08
			
			// number of lines read
			line.get(input_lines, false, size_t(5e05)); // default: 0.5 million lines
			
			// locally save tmp files
			line.get(local_tmp_files, false, false);
			
			// ancestral alleles
			line.get(input_anc_file, false);
			
			// remove unmatched ancestral alleles
			line.get(require_anc_match, false, false);
			
			// require SNPs
			line.get(require_snp, false, true); // default: only SNPs
		}
		else
		{
			// target position
			const bool do_position = line.get(share_position, false);
			
			// batch file of positions
			const bool do_batch = line.get(share_batch, false);
			
			if (do_position && do_batch)
				throw std::invalid_argument("Conflicting target position input");
			
			if (!do_position && !do_batch)
				throw std::invalid_argument("Missing target position input");
			
			
			line.get(input_hmm_file, true);
			
			line.get(effective_size, false, size_t(10000)); // default: 10000
			line.get(mutation_rate, false, double(1e-08)); // default: 1e-08
			line.get(max_missing, false, double(0.05)); // defaul: 5%
			line.get(age_limit_sharers, false, size_t(100)); // default: 100
			line.get(age_outgroup_size, false, size_t(100)); // default: 100
		}

		line.finish();
	}
	catch(const int & help)
	{
		return EXIT_SUCCESS;
	}
	catch(const std::exception & error)
	{
		std::cout << error.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	
	// set seed
	if (seed.good())
	{
		set_random_seed(seed);
	}
	
	
	// redirect log/err output
	Redirect redirect_log(std::clog, output.value + ".log");
	Redirect redirect_err(std::cerr, output.value + ".err");
	
	
	
	std::cout << std::endl;
	std::cout << "Genealogical estimation of variant age (GEVA)" << std::endl;
	std::cout << "(c) Patrick K. Albers" << std::endl;
	std::cout << std::endl;
	
	
	// print programme call to log
	line.print(std::clog);
	
	
	// begin timer
	Clock runtime;
	
	
	// main algorithm
	try
	{
		// process input data
		
		Gen::Grid::Data grid;
		
		if (!input_bin_file.good())
		{
			// switch between map rec rate or fixed rec rate
			Gen::Map gmap = (input_map_file.good()) ? load_map(input_map_file): load_map(input_rec_rate);
			
			// ancestral alleles
			std::string amap;
			if (input_anc_file.good())
				amap = input_anc_file;
			
			input_bin_file = load_vcf(input_vcf_file,
																gmap,
																amap,
																output,
																input_lines,
																false,
																local_tmp_files,
																require_anc_match,
																require_snp);
			
			grid = load_bin(input_bin_file);
			
			grid->print_sample(output.value + ".sample.txt");
			grid->print_marker(output.value + ".marker.txt");
		}
		else
		{
			grid = load_bin(input_bin_file);
			grid->cache(buffer);
			
			Gen::Share::Data share = std::make_shared<Gen::Share>();
			
			if (share_position.good())
				select_share(share, grid, share_position);
			
			if (share_batch.good())
				select_share(share, grid, share_batch);
			
			
			// prepare hmm
			
			IBD::DetectMethod method = IBD::DETECT_HMM;
			IBD::HMM::Model::Data hmm_model;
			hmm_model = load_hmm(grid, input_hmm_file[0], input_hmm_file[1], effective_size, output);
			
			
			// execute age estimation
			
			Age::Param::Data param = std::make_shared< Age::Param >(grid, effective_size, mutation_rate);
			
			if (age_limit_sharers.good())
				param->limit_sharers = age_limit_sharers.value;
			
			if (age_outgroup_size.good())
				param->outgroup_size = age_outgroup_size.value;
			
			param->threads = thread;
			
			infer_age(param, method, max_missing, output, false, false, share, grid, runtime, hmm_model, nullptr, thread);
		}
	}
	catch(const std::exception & error)
	{
		std::cout << error.what() << std::endl;
		return EXIT_FAILURE;
	}
	
	
	std::cout << "This is the end!" << std::endl;
	std::clog << "This is the end!" << std::endl;
	runtime.print(std::cout);
	runtime.print(std::clog);
	
	
	return EXIT_SUCCESS;
}

