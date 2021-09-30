//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef load_vcf_h
#define load_vcf_h


#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Clock.hpp"
#include "Progress.hpp"
#include "LoadVcf.hpp"
#include "Gen.hpp"
#include "GenGrid.hpp"



static void make_anc_allele_map(const std::string & filename, LoadVcf::AncAlleleMap & amap)
{
	std::ifstream file(filename);
	std::string   line;

	while (std::getline(file, line))
	{
		std::string chr_str;
		char *      chr_ptr = nullptr;
		int         chr = -1;
		size_t      pos = 0;
		std::string anc;

		std::istringstream iss(line);

		if (!(iss >> chr_str >> pos >> anc))
		{
			file.close();
			throw std::runtime_error(std::string("Unable to read file:  ") + filename);
		}

		if (chr_str.size() > 3 &&
				(chr_str[0] == 'c' || chr_str[0] == 'C') &&
				(chr_str[1] == 'h' || chr_str[1] == 'H')
				(chr_str[2] == 'r' || chr_str[2] == 'R'))
		{
			if (chr_str[3] == 'X')
			{
				chr = 23;
			}
			else
			{
				chr_ptr = &chr_str[3];
			}
		}
		else
		{
			if (chr_str.size() == 1 && chr_str[0] == 'X')
			{
				chr = 23;
			}
			else
			{
				chr_ptr = &chr_str[0];
			}
		}

		if (chr_ptr != nullptr)
		{
			char * err;
			float  val = std::strtof(chr_ptr, &err);
			if (*err != '\0')
				throw std::invalid_argument("Cannot parse chromosome identifier");
			chr = static_cast<int>(val);
		}

		amap[chr][pos] = anc;
	}

	file.close();
}


inline std::string load_vcf(const std::string & input,
														const Gen::Map & gmap,
														const std::string & amap,
														const std::string & outfile,
														const size_t & buffer_limit,
														const bool compress = false,
														const bool local_tmp = false,
														const bool require_anc_match = false,
														const bool require_snp = true)
//							const size_t chunk_beg = 0, const size_t chunk_end = 0)
{
	std::cout << "Loading variant data from VCF file to binary file" << std::endl << "<< " << input << std::endl;
	std::clog << "Loading variant data from VCF file to binary file" << std::endl << "<< " << input << std::endl;


	const std::string grid_file = outfile + ".bin"; // filename

	std::cout << ">> " << grid_file << std::endl;
	std::clog << ">> " << grid_file << std::endl;

	std::cout << std::endl;
	std::clog << std::endl;


	size_t warn = 0; // count warnings

	try
	{
		Clock   time;
		LoadVcf load(input); // reads sample vector

		if (!amap.empty())
		{
			make_anc_allele_map(amap, load.ancestral);
			load.has_ancestral = true;
			load.filter.remove_unmatched_ancestral = require_anc_match;

			std::cout << " Determine ancestral allelic states according to file:  " << amap << std::endl;
			std::clog << " Determine ancestral allelic states according to file:  " << amap << std::endl;

			if (load.filter.remove_unmatched_ancestral)
			{
				std::cout << " Removing variants where the ancestral allele does not match given ref/alt alleles" << std::endl;
				std::clog << " Removing variants where the ancestral allele does not match given ref/alt alleles" << std::endl;
			}
		}
		else
			load.has_ancestral = false;

		load.filter.require_snp = require_snp;
		if (load.filter.require_snp)
		{
			std::cout << " Retaining only SNPs" << std::endl;
			std::clog << " Retaining only SNPs" << std::endl;
		}

		const size_t sample_size = load.sample.size();
		const size_t buffer_size = static_cast<size_t>(static_cast<double>(sample_size * buffer_limit) / static_cast<double>(1024 * 1024 * sizeof(Gen::value_t)));

		std::cout << " Detected sample size: " << sample_size << " individuals" << std::endl;
		std::clog << " Detected sample size: " << sample_size << " individuals" << std::endl;

		std::cout << " Buffer window size: " << buffer_limit << " variants (~" << buffer_size << " Mb)" << std::endl;
		std::clog << " Buffer window size: " << buffer_limit << " variants (~" << buffer_size << " Mb)" << std::endl;

//		std::cout << " Data compression: " << (compress ? "On": "Off") << std::endl;
//		std::clog << " Data compression: " << (compress ? "On": "Off") << std::endl;

//		if (chunk_beg != chunk_end)
//		{
//			load.filter.position_beg = std::min(chunk_beg, chunk_end);
//			load.filter.position_end = std::max(chunk_beg, chunk_end);
//
//			std::cout << " Chunk: positions from " << load.filter.position_beg << " (inclusive) to " << load.filter.position_end << " (exclusive)" << std::endl;
//			std::clog << " Chunk: positions from " << load.filter.position_beg << " (inclusive) to " << load.filter.position_end << " (exclusive)" << std::endl;
//		}

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


#endif /* load_vcf_h */
