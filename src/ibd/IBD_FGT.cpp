//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "IBD_FGT.hpp"


using namespace Gen;
using namespace IBD;
using namespace FGT;


// Gamete frequencies

// construct

Frequency::Frequency(const Marker::Key & focal, const Sample::Key::Vector & sharers, const Grid::Data source)
: grid(source)
, freq(source->marker_size(), 0)
{
	size_t k = 0;
	
	Sample::Key::Vector::const_iterator it, ti = sharers.cend();
	
	for (it = sharers.cbegin(); it !=ti; ++it)
	{
		const Variant::Vector::Data var = this->grid->get(*it);
		
		const ChrType chr = identify(var, focal);
		
		if (chr == CHR_VOID)
		{
			throw std::logic_error("Unable to identify shared haplotype at focal site");
		}
		
		for (size_t i = 0, n = this->grid->marker_size(); i < n; ++i)
		{
			const hap_t hap = var->hap(chr, i);
			
			if (is_haplotype<H_>(hap))
			{
				this->freq[i] = this->grid->sample_size() * 2; // to skip missing/unknown sites
				continue;
			}
			
			if (is_haplotype<H1>(hap))
			{
				this->freq[i] += 1;
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


// identify focal haplotype

ChrType Frequency::identify(const Variant::Vector::Data var, const Marker::Key & site)
{
	Variant foc = var->at(site);
	
	const HapType mat = haplotype_to_index(foc.hap(MATERNAL));
	const HapType pat = haplotype_to_index(foc.hap(PATERNAL));
	
	if (is_haplotype<H1>(mat) && is_haplotype<H0>(pat))
	{
		return MATERNAL;
	}
	
	if (is_haplotype<H1>(pat) && is_haplotype<H0>(mat))
	{
		return PATERNAL;
	}
	
	return CHR_VOID;
}





// Four-Gamates algorithm

// construct

Algorithm::Algorithm(const Variant::Vector::Data v0, const Variant::Vector::Data v1) // classic
: sharers(v0, v1)
, chromos(CHR_VOID, CHR_VOID)
, gametes(nullptr)
, length(v0->size())
{
	if (this->length == 0)
		throw std::runtime_error("Variant vectors are of zero size");
	
	if (this->length != v1->size())
		throw std::runtime_error("Variant vectors are of unequal size");
	
	if (!v0->is_phased() || !v1->is_phased())
		throw std::runtime_error("Variant vector is not phased");
}

Algorithm::Algorithm(const Variant::Vector::Data v0, const ChrType c0, const Variant::Vector::Data v1, const ChrType c1, const Frequency::Data freq) // frequency dependent
: sharers(v0, v1)
, chromos(c0, c1)
, gametes(freq)
, length(v0->size())
{
	if (this->length == 0)
		throw std::runtime_error("Variant vectors are of zero size");
	
	if (this->length != v1->size())
		throw std::runtime_error("Variant vectors are of unequal size");
	
	if (!v0->is_phased() || !v1->is_phased())
		throw std::runtime_error("Variant vector is not phased");
}


// detect segment

Segment Algorithm::detect(const Marker::Key & focal)
{
	if (this->gametes)
	{
		return this->detect_frqbase(focal);
	}
	return this->detect_classic(focal);
}


// detection methods

inline Algorithm::four_gam_t four_gametes(const Variant & a, const Variant & b)
{
	return {{
		haplotype_to_index(a.hap(MATERNAL)),
		haplotype_to_index(a.hap(PATERNAL)),
		haplotype_to_index(b.hap(MATERNAL)),
		haplotype_to_index(b.hap(PATERNAL))
	}};
}


inline bool four_gametes_test(const Algorithm::four_gam_t & focal, const Algorithm::four_gam_t & trans)
{
	bool test[2][2] = { {false, false}, {false, false} };
	
	for (unsigned i = 0; i < Algorithm::four_gam_n; ++i)
	{
		if (is_haplotype<H_>(focal[i])) return false;
		if (is_haplotype<H_>(trans[i])) return false;
		
		const unsigned f = static_cast<unsigned>(is_haplotype<H1>(focal[i]));
		const unsigned t = static_cast<unsigned>(is_haplotype<H1>(trans[i]));

		test[f][t] = true;
	}
	
	return (test[0][0] && test[0][1] && test[1][0] && test[1][1]);
}

Segment Algorithm::detect_classic(const Marker::Key & focal)
{
	Variant foc0 = this->sharers.first->at(focal);
	Variant foc1 = this->sharers.second->at(focal);
	
	const four_gam_t focal_four = four_gametes(foc0, foc1);
	
	
	Segment segment(0, this->length - 1); // initialise to full range
	
	
	// LHS
	
	for (size_t i = focal.value; i > 0; --i)
	{
		const size_t j = i - 1;
		
		Variant trans0 = this->sharers.first->at(j);
		Variant trans1 = this->sharers.second->at(j);
		
		const four_gam_t trans_four = four_gametes(trans0, trans1);
		
		if (four_gametes_test(focal_four, trans_four))
		{
			segment[LHS] = j;
			break;
		}
	}
	
	
	// RHS
	
	for (size_t i = focal.value; i < this->length - 1; ++i)
	{
		const size_t j = i + 1;
		
		Variant trans0 = this->sharers.first->at(j);
		Variant trans1 = this->sharers.second->at(j);
		
		const four_gam_t trans_four = four_gametes(trans0, trans1);
		
		if (four_gametes_test(focal_four, trans_four))
		{
			segment[RHS] = j;
			break;
		}
	}
	
	return segment;
}


Segment Algorithm::detect_frqbase(const Marker::Key & focal)
{
	const size_t focal_frq = this->gametes->sample(focal);
	
	Segment segment(0, this->length - 1); // initialise to full range
	
	
	// LHS
	
	for (size_t i = focal.value; i > 0; --i)
	{
		const size_t j = i - 1;
		
		const hap_t hap0 = this->sharers.first->at(j).hap( this->chromos.first );
		const hap_t hap1 = this->sharers.second->at(j).hap( this->chromos.second );
		
		if (hap_break(hap0, hap1) == HAP_BREAK)
		{
			const size_t sample_frq = this->gametes->sample(j);
			
			if (sample_frq == 1)
			{
				continue;
			}
			
			const size_t subset_frq = this->gametes->subset(j);
			
			if (subset_frq < focal_frq && subset_frq < sample_frq)
			{
				segment[LHS] = j;
				break;
			}
		}
	}
	
	
	// RHS
	
	for (size_t i = focal.value; i < this->length - 1; ++i)
	{
		const size_t j = i + 1;
		
		const hap_t hap0 = this->sharers.first->at(j).hap( this->chromos.first );
		const hap_t hap1 = this->sharers.second->at(j).hap( this->chromos.second );
		
		if (hap_break(hap0, hap1) == HAP_BREAK)
		{
			const decimal_t sample_frq = this->gametes->sample(j);
			
			if (sample_frq == 1)
			{
				continue;
			}
			
			const decimal_t subset_frq = this->gametes->subset(j);
			
			if (subset_frq < focal_frq && subset_frq < sample_frq)
			{
				segment[RHS] = j;
				break;
			}
		}
	}
	
	
	return segment;
}


// print to stream

void Algorithm::print(const Gen::Sample::Key::Pair & pair, const Gen::Marker::Key & focal, std::ostream & stream) const
{
	Variant foc0 = this->sharers.first->at(focal);
	Variant foc1 = this->sharers.second->at(focal);
	
	const four_gam_t focal_four = four_gametes(foc0, foc1);
	
	
	// LHS
	
	for (size_t i = focal.value; i > 0; --i)
	{
		const size_t j = i - 1;
		
		Variant trans0 = this->sharers.first->at(j);
		Variant trans1 = this->sharers.second->at(j);
		
		const four_gam_t trans_four = four_gametes(trans0, trans1);
		
		if (four_gametes_test(focal_four, trans_four))
		{
			stream << pair.first << ' ' << pair.second << ' ' << focal.value << ' ' << j << std::endl;
		}
	}
	
	
	// RHS
	
	for (size_t i = focal.value; i < this->length - 1; ++i)
	{
		const size_t j = i + 1;
		
		Variant trans0 = this->sharers.first->at(j);
		Variant trans1 = this->sharers.second->at(j);
		
		const four_gam_t trans_four = four_gametes(trans0, trans1);
		
		if (four_gametes_test(focal_four, trans_four))
		{
			stream << pair.first << ' ' << pair.second << ' ' << focal.value << ' ' << j << std::endl;
		}
	}
}

void Algorithm::print_header(std::ostream & stream)
{
	stream << "SampleID0 SampleID1 FocalID BreakID" << std::endl;
}

