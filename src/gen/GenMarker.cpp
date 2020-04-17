//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "GenMarker.hpp"

#include <cstring>


using namespace Gen;


// constructs

Marker::Marker()
: index(0)
, chromosome(-1)
, position(0)
, hap_count() // default initialise !!!
, gen_count() // default initialise !!!
, rec_rate(-1)
, gen_dist(-1)
{}

/*
Marker::Marker(const Marker & other)
: index(other.index)
, label(other.label)
, chromosome(other.chromosome)
, position(other.position)
, allele(other.allele)
, hap_count(other.hap_count)
, gen_count(other.gen_count)
, rec_rate(other.rec_rate)
, gen_dist(other.gen_dist)
{}

Marker::Marker(Marker && other)
: index(other.index)
, label(std::move(other.label))
, chromosome(other.chromosome)
, position(other.position)
, allele(other.allele)
, hap_count(std::move(other.hap_count))
, gen_count(std::move(other.gen_count))
, rec_rate(other.rec_rate)
, gen_dist(other.gen_dist)
{}
*/


// count alleles and genotype

void Marker::count(const gen_t g)
{
	const hap_pair_t hh = genotype_to_haplotypes(g);
	
	const HapType mat = haplotype_to_index(hh[ MATERNAL ]);
	const HapType pat = haplotype_to_index(hh[ PATERNAL ]);
	
	switch (mat)
	{
		case H0:
		{
			++this->hap_count[H0];
			
			switch (pat)
			{
				case H0: ++this->hap_count[H0]; ++this->gen_count[G0]; break;
				case H1: ++this->hap_count[H1]; ++this->gen_count[G1]; break;
				case H_: ++this->hap_count[H_]; ++this->gen_count[G_]; break;
			}
			break;
		}
		case H1:
		{
			++this->hap_count[H1];
			
			switch (pat)
			{
				case H0: ++this->hap_count[H0]; ++this->gen_count[G1]; break;
				case H1: ++this->hap_count[H1]; ++this->gen_count[G2]; break;
				case H_: ++this->hap_count[H_]; ++this->gen_count[G_]; break;
			}
			break;
		}
		case H_:
		{
			++this->hap_count[H_];
			
			switch (pat)
			{
				case H0: ++this->hap_count[H0]; ++this->gen_count[G_]; break;
				case H1: ++this->hap_count[H1]; ++this->gen_count[G_]; break;
				case H_: ++this->hap_count[H_]; ++this->gen_count[G_]; break;
			}
			break;
		}
	}
}


// print to stream

void Marker::print(std::ostream & stream) const
{
	stream << this->index.value << ' ';
	stream << this->chromosome << ' ';
	stream << this->position << ' ';
	stream << this->label << ' ';
	stream << this->allele.str() << ' ';
	stream << this->hap_count[H0] << ' '  << this->hap_count[H1] << ' ' << this->hap_count[H_] << ' ';
	stream << this->gen_count[G0] << ' '  << this->gen_count[G1] << ' ' << this->gen_count[G2] << ' ' << this->gen_count[G_] << ' ';
	stream << std::fixed << std::setprecision(8) << this->rec_rate << ' ';
	stream << std::fixed << std::setprecision(8) << this->gen_dist;
}


// return as string

std::string Marker::str() const
{
	std::ostringstream oss;
	this->print(oss);
	return oss.str();
}


// header for printing

const std::string Marker::header = "MarkerID Chromosome Position Label Alleles "
                                   "AlleleCount0 AlleleCount1 AlleleCountX "
                                   "GenotypeCount0 GenotypeCount1 GenotypeCount2 GenotypeCountX "
                                   "RecombRate GenDist";



// AlleleList

// constructs

Marker::Allele::Allele()
: value(-1)
{}

/*
Marker::Allele::Allele(const Allele & other)
: value(other.value)
, cache(other.cache)
{}

Marker::Allele::Allele(Allele && other)
: value(other.value)
, cache(std::move(other.cache))
{}
*/


// parse allele string
void Marker::Allele::parse(const std::string & raw)
{
	if (raw.empty())
	{
		return;
	}
	
	// determine SNP
	if (raw.size() == 3 && raw[1] == Marker::Allele::sep)
	{
		this->value = Marker::Allele::mask_snp(raw[0], raw[2]);
	}
	
	// store formated string
	if (this->value < 0)
	{
		const char * p = raw.data();
		const char * q = std::strchr(p, ',');
		while (q)
		{
			if (q-p > 0)
				this->cache.push_back(std::string(p, q-p));
			
			p = q + 1;
			q = std::strchr(p, ',');
		}
		this->cache.push_back(std::string(p));
		
//		std::vector<char> chunk(raw.size());
//
//		chunk.reserve(256);
//
//		std::string::const_iterator it, end = raw.cend();
//
//		for (it = raw.cbegin(); it != end; ++it)
//		{
//			const char c = *it;
//
//			if (isalpha(c))
//			{
//				chunk.push_back(*it);
//				continue;
//			}
//
//			if (*it == Marker::Allele::sep || it + 1 == end)
//			{
//				if (!chunk.empty())
//				{
//					this->cache.push_back( std::string(chunk.begin(), chunk.end()) );
//					chunk.clear();
//				}
//			}
//		}
	}
	
	// re-evaluate if SNP
//	if (this->cache.size() == 2 && this->cache[0].size() == 1 && this->cache[1].size() == 1)
//	{
//		this->value = Marker::Allele::mask_snp(this->cache[0][0], this->cache[1][0]);
//
//		if (this->value >= 0)
//		{
//			this->cache = cache_t(0);
//		}
//	}
}


// return alleles

std::string Marker::Allele::operator [] (const size_t i) const
{
	static const std::array< char, 2 > mask[25] = { // output mask to return usual alleles
		{'.','.'},  {'.','A'},  {'.','C'},  {'.','G'},  {'.','T'},
		{'A','.'},  {'A','A'},  {'A','C'},  {'A','G'},  {'A','T'},
		{'C','.'},  {'C','A'},  {'C','C'},  {'C','G'},  {'C','T'},
		{'G','.'},  {'G','A'},  {'G','C'},  {'G','G'},  {'G','T'},
		{'T','.'},  {'T','A'},  {'T','C'},  {'T','G'},  {'T','T'}
	};
	
	if (this->value < 0)
	{
		return this->cache.at(i);
	}
	
	return std::string(1, mask[ this->value ].at(i));
}


// check if bi-allelic SNP

bool Marker::Allele::is_snp() const
{
	return (this->value >= 0);
}


// return number of alleles

size_t Marker::Allele::size() const
{
	if (this->value < 0)
	{
		return this->cache.size();
	}
	
	return 2;
}


// return raw string

std::string Marker::Allele::str() const
{
	std::ostringstream oss;
	
	if (this->value < 0)
	{
		const size_t size = this->cache.size();
		
		if (size == 0)
		{
			return "???";
		}
		
		for (size_t i = 0; i < size; ++i)
		{
			oss << this->cache[i];
			
			if (i != size - 1)
			{
				oss << Marker::Allele::sep;
			}
		}
	}
	else
	{
		oss << this->operator[](0);
		oss << Marker::Allele::sep;
		oss << this->operator[](1);
	}
	
	return oss.str();
}


// encode SNP to determine mask value

Marker::Allele::value_t Marker::Allele::mask_snp(char c0, char c1)
{
	c0 = toupper(c0);
	c1 = toupper(c1);
	
	if (isalpha(c0) && isalpha(c1))
	{
		switch (c0)
		{
			case '.': switch (c1) {
				case '.': return  0; break;
				case 'A': return  1; break;
				case 'C': return  2; break;
				case 'G': return  3; break;
				case 'T': return  4; break;
				default: break; } break;
			case 'A': switch (c1) {
				case '.': return  5; break;
				case 'A': return  6; break;
				case 'C': return  7; break;
				case 'G': return  8; break;
				case 'T': return  9; break;
				default: break; } break;
			case 'C': switch (c1) {
				case '.': return 10; break;
				case 'A': return 11; break;
				case 'C': return 12; break;
				case 'G': return 13; break;
				case 'T': return 14; break;
				default: break; } break;
			case 'G': switch (c1) {
				case '.': return 15; break;
				case 'A': return 16; break;
				case 'C': return 17; break;
				case 'G': return 18; break;
				case 'T': return 19; break;
				default: break; } break;
			case 'T': switch (c1) {
				case '.': return 20; break;
				case 'A': return 21; break;
				case 'C': return 22; break;
				case 'G': return 23; break;
				case 'T': return 24; break;
				default: break; } break;
			default: break;
		}
	}
	
	return -1;
}

