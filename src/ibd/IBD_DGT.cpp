//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "IBD_DGT.hpp"


using namespace Gen;
using namespace IBD;
using namespace DGT;



// Allele frequencies

// construct

Frequency::Frequency(const Sample::Key::Vector & sharers, const Grid::Data source)
: grid(source)
, freq(source->marker_size(), 0)
{
	size_t k = 0;
	
	Sample::Key::Vector::const_iterator it, ti = sharers.cend();
	
	for (it = sharers.cbegin(); it !=ti; ++it)
	{
		const Variant::Vector::Data var = this->grid->get(*it);
		
		for (size_t i = 0, n = this->grid->marker_size(); i < n; ++i)
		{
			const gen_t gen = var->gen(i);
			
			if (is_genotype<G_>(gen))
			{
				this->freq[i] = this->grid->sample_size() * 2; // to skip missing/unknown sites
				continue;
			}
			
			if (is_genotype<G1>(gen))
			{
				this->freq[i] += 1;
				continue;
			}
			
			if (is_genotype<G2>(gen))
			{
				this->freq[i] += 2;
			}
		}
		
		++k;
	}
	
	if (k == 0)
	{
		throw std::logic_error("Empty set of sharers");
	}
}


// return frequencies

size_t Frequency::sample(const Marker::Key & site) const
{
	return this->grid->marker(site).hap_count[H1];
}

size_t Frequency::subset(const Marker::Key & site) const
{
	return this->freq.at(site.value);
}



// construct

Algorithm::Algorithm(const Gen::Variant::Vector::Data v0, const Gen::Variant::Vector::Data v1)
: sharers(v0, v1)
, alleles(nullptr)
, length(v0->size())
{
	if (this->length == 0)
		throw std::runtime_error("Variant vectors are of zero size");
	
	if (this->length != v1->size())
		throw std::runtime_error("Variant vectors are of unequal size");
	
	
//	if (a->size() != b->size())
//	{
//		throw std::runtime_error("Variant vectors are of unequal size");
//	}
//	
//	gen_vector_t const & g0 = a->gen();
//	gen_vector_t const & g1 = b->gen();
//
//	
//	// find discordances
//	
//	for (size_t i = 0; i < this->length; ++i)
//	{
//		const GenType i0 = genotype_to_index(g0[i]);
//		const GenType i1 = genotype_to_index(g1[i]);
//		
//		if ((is_genotype<G0>(i0) && is_genotype<G2>(i1)) ||
//			(is_genotype<G2>(i0) && is_genotype<G0>(i1)))
//		{
//			this->brk.insert(i);
//		}
//	}
//	
//	// insert bounds
//	this->brk.insert(0);
//	this->brk.insert(this->length - 1);
}

Algorithm::Algorithm(const Variant::Vector::Data v0, const Variant::Vector::Data v1, const Frequency::Data freq) // frequency dependent
: sharers(v0, v1)
, alleles(freq)
, length(v0->size())
{
	if (this->length == 0)
		throw std::runtime_error("Variant vectors are of zero size");
	
	if (this->length != v1->size())
		throw std::runtime_error("Variant vectors are of unequal size");
}


// detect segment

Segment Algorithm::detect(const Marker::Key & focal)
{
	if (this->alleles)
	{
		return this->detect_frqbase(focal);
	}
	return this->detect_classic(focal);
}

Segment Algorithm::detect_classic(const Marker::Key & focal)
{
	Segment segment(0, this->length - 1); // initialise to full range
	
	
	// LHS
	
	for (size_t i = focal.value; i > 0; --i)
	{
		const size_t j = i - 1;
		
		const gen_t gen0 = this->sharers.first->at(j).gen();
		const gen_t gen1 = this->sharers.second->at(j).gen();
		
		if ((is_genotype<G0>(gen0) && is_genotype<G2>(gen1)) ||
			(is_genotype<G2>(gen0) && is_genotype<G0>(gen1)))
		{
			segment[LHS] = j;
			break;
		}
	}
	
	
	// RHS
	
	for (size_t i = focal.value; i < this->length - 1; ++i)
	{
		const size_t j = i + 1;
		
		const gen_t gen0 = this->sharers.first->at(j).gen();
		const gen_t gen1 = this->sharers.second->at(j).gen();
		
		if ((is_genotype<G0>(gen0) && is_genotype<G2>(gen1)) ||
			(is_genotype<G2>(gen0) && is_genotype<G0>(gen1)))
		{
			segment[RHS] = j;
			break;
		}
	}
	
	return segment;
	
//	discord_set_t::const_iterator upper;
//	discord_set_t::const_iterator lower;
//	
//	if (site.value > 0)
//	{
//		upper = this->brk.lower_bound(site.value);
//	}
//	else
//	{
//		upper = this->brk.upper_bound(site.value);
//	}
//	
//	lower = std::prev(upper);
//	
//	return Segment(*lower, *upper);
}

Segment Algorithm::detect_frqbase(const Marker::Key & focal)
{
	const size_t focal_frq = this->alleles->sample(focal);
	
	Segment segment(0, this->length - 1); // initialise to full range
	
	
	// LHS
	
	for (size_t i = focal.value; i > 0; --i)
	{
		const size_t j = i - 1;
		
		const gen_t gen0 = this->sharers.first->at(j).gen();
		const gen_t gen1 = this->sharers.second->at(j).gen();
		
		if (gen_break(gen0, gen1) == GEN_BREAK)
		{
			const size_t sample_frq = this->alleles->sample(j);
			
//			if (sample_frq == 1)
//			{
//				continue;
//			}
			
			const size_t subset_frq = this->alleles->subset(j);
			
			if (subset_frq < focal_frq && subset_frq < sample_frq)
			{
				segment[LHS] = i;
				break;
			}
		}
	}
	
	
	// RHS
	
	for (size_t i = focal.value; i < this->length - 1; ++i)
	{
		const size_t j = i + 1;
		
		const gen_t gen0 = this->sharers.first->at(j).gen();
		const gen_t gen1 = this->sharers.second->at(j).gen();
		
		if (gen_break(gen0, gen1) == GEN_BREAK)
		{
			const decimal_t sample_frq = this->alleles->sample(j);
			
//			if (sample_frq == 1)
//			{
//				continue;
//			}
			
			const decimal_t subset_frq = this->alleles->subset(j);
			
			if (subset_frq < focal_frq && subset_frq < sample_frq)
			{
				segment[RHS] = i;
				break;
			}
		}
	}
	
	
	return segment;
}



// print to stream

void Algorithm::print(const Gen::Sample::Key::Pair & pair, const Gen::Marker::Key & focal, std::ostream & stream) const
{
//	const size_t last = *this->brk.crbegin();
//	
//	discord_set_t::const_iterator it, ti = this->brk.cend();
//	
//	for (it = this->brk.cbegin(); it != ti; ++it)
//	{
//		const size_t x = *it;
//		
//		if (x == 0 || x == last)
//			continue;
//		
//		stream << pair.first << ' ' << pair.second << ' ' << site.value << ' ' << (x - 1) << std::endl;
//	}

	
	// LHS
	
	for (size_t i = focal.value; i > 0; --i)
	{
		const size_t j = i - 1;
		
		const gen_t gen0 = this->sharers.first->at(j).gen();
		const gen_t gen1 = this->sharers.second->at(j).gen();
		
		if ((is_genotype<G0>(gen0) && is_genotype<G2>(gen1)) ||
			(is_genotype<G2>(gen0) && is_genotype<G0>(gen1)))
		{
			stream << pair.first << ' ' << pair.second << ' ' << focal.value << ' ' << j << std::endl;
		}
	}
	
	
	// RHS
	
	for (size_t i = focal.value; i < this->length - 1; ++i)
	{
		const size_t j = i + 1;
		
		const gen_t gen0 = this->sharers.first->at(j).gen();
		const gen_t gen1 = this->sharers.second->at(j).gen();
		
		if ((is_genotype<G0>(gen0) && is_genotype<G2>(gen1)) ||
			(is_genotype<G2>(gen0) && is_genotype<G0>(gen1)))
		{
			stream << pair.first << ' ' << pair.second << ' ' << focal.value << ' ' << j << std::endl;
		}
	}
}

void Algorithm::print_header(std::ostream & stream)
{
	stream << "SampleID0 SampleID1 FocalID BreakID" << std::endl;
}
