//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBD_hpp
#define IBD_hpp

#include <memory>
#include <utility>

#include "Decimal.h"

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenVariant.hpp"
#include "GenGrid.hpp"


namespace IBD
{
	// Detection methods
	enum DetectMethod : unsigned
	{
		DETECT_DGT = 0,
		DETECT_FGT = 1,
		DETECT_HMM = 2,
		DETECT_SIM = 3,
		DETECT_VOID = 4 // undefined
	};
	
	
	// Left/right hand side
	enum LR : uint_fast8_t
	{
		LHS = 0,
		RHS = 1
	};
	
	
	// Generic left/right hand side container
	template<typename T>
	class Side
	{
	public:
		
		using Type = T;
		
		// constructs
		Side() : lhs(), rhs() {}
		Side(const T & l, const T & r) : lhs(l) , rhs(r) {}
		Side(T && l, T && r) : lhs(std::move(l)) , rhs(std::move(r)) {}
		Side(const Side<T> & other) : lhs(other.lhs) , rhs(other.rhs) {}
		Side(Side<T> && other) : lhs(std::move(other.lhs)) , rhs(std::move(other.rhs)) {}
		
		// assigns
		Side<T> & operator = (const Side<T> & other) { this->lhs = other.lhs; this->rhs = other.rhs; return *this; }
		Side<T> & operator = (Side<T> && other) { this->lhs = std::move(other.lhs); this->rhs = std::move(other.rhs); return *this; }
		
		// refer
		constexpr T const & operator [] (const LR s) const { return (s == LHS) ? this->lhs: this->rhs; }
		T & operator [] (const LR s) { return (s == LHS) ? this->lhs: this->rhs; }
		
		// sort/compare other segment
		bool operator <  (const Side<T> & other) const { return (this->lhs == other.lhs) ? (this->rhs < other.rhs): (this->lhs < other.lhs); }
		bool operator  > (const Side<T> & other) const { return (this->lhs == other.lhs) ? (this->rhs > other.rhs): (this->lhs > other.lhs); }
		bool operator == (const Side<T> & other) const { return (this->lhs == other.lhs && this->rhs == other.rhs); }
		bool operator != (const Side<T> & other) const { return (this->lhs != other.lhs || this->rhs != other.rhs); }
		
		
	private:
		
		T lhs;
		T rhs;
	};
	
	// IBD segment, inclusive endpoints on each side
	using Segment = Side< Gen::Marker::Key >;
	
	// Mutational differences over haplotypes per side
	using SegDiff = Side< int >;
	
	// Distributions on each side
	using Distribution = Side< decimal_vector_t >;
	
	
	inline Distribution init_distribution(const size_t & n0, const size_t & n1, const decimal_t value = decimal_nil)
	{
		return Distribution(decimal_vector_t(n0, value), decimal_vector_t(n1, value));
	}
	
	
	
	// Calculate missing rate in pairwise runs
	
	inline decimal_t missing_rate(const Gen::gen_vector_t & a, const Gen::gen_vector_t & b, const size_t & n)
	{
		size_t r = 0;
		
		for (size_t i = 0; i < n; ++i)
		{
			if (Gen::is_genotype<Gen::G_>(a.at(i)) || Gen::is_genotype<Gen::G_>(b.at(i))) ++r;
		}
		
		return static_cast<decimal_t>(r) / static_cast<decimal_t>(n);
	}
	
	
	
	// Haplotype breakpoints
	
	enum HapBreak : uint_fast8_t
	{
		HAP_SHARE = 0,
		HAP_BREAK = 1,
		HAP_UNDEF = 2  // undefined
	};
	
	using brk_vector_t = std::vector< HapBreak >;
	
	
	inline HapBreak hap_break(const Gen::hap_t h0, const Gen::hap_t h1)
	{
		using namespace Gen;
		
		const HapType i0 = haplotype_to_index(h0);
		const HapType i1 = haplotype_to_index(h1);
		
		if (is_haplotype<H_>(i0) || is_haplotype<H_>(i1))
		{
			return HAP_UNDEF;
		}
		
		return (i0 == i1) ? HAP_SHARE: HAP_BREAK;
	}


	// Genotype breakpoints
	
	enum GenBreak : uint_fast8_t
	{
		GEN_SHARE = 0,
		GEN_BREAK = 1,
		GEN_UNDEF = 2  // undefined
	};
	
	using brk_vector_t = std::vector< HapBreak >;
	
	
	inline GenBreak gen_break(const Gen::gen_t g0, const Gen::gen_t g1)
	{
		using namespace Gen;
		
		const GenType i0 = genotype_to_index(g0);
		const GenType i1 = genotype_to_index(g1);
		
		if (is_genotype<G0>(i0))
		{
			if (is_genotype<G0>(i1)) return GEN_SHARE;
			if (is_genotype<G1>(i1)) return GEN_BREAK;
			if (is_genotype<G2>(i1)) return GEN_BREAK;
		}
		if (is_genotype<G1>(i0))
		{
			if (is_genotype<G0>(i1)) return GEN_BREAK;
			if (is_genotype<G1>(i1)) return GEN_UNDEF;
			if (is_genotype<G2>(i1)) return GEN_BREAK;
		}
		if (is_genotype<G2>(i0))
		{
			if (is_genotype<G0>(i1)) return GEN_BREAK;
			if (is_genotype<G1>(i1)) return GEN_BREAK;
			if (is_genotype<G2>(i1)) return GEN_SHARE;
		}
		
		return GEN_UNDEF;
	}
	
	
	inline std::pair<Gen::ChrType, Gen::ChrType> select_shared_chromosomes(const Gen::gen_t g0, const Gen::gen_t g1)
	{
		std::pair<Gen::ChrType, Gen::ChrType> chr;
		
		chr.first  = Gen::CHR_VOID;
		chr.second = Gen::CHR_VOID;
		
		Gen::hap_pair_t hh0 = Gen::genotype_to_haplotypes(g0);
		Gen::hap_pair_t hh1 = Gen::genotype_to_haplotypes(g1);
		
		if (Gen::is_haplotype<Gen::H1>(hh0[0]) && Gen::is_haplotype<Gen::H0>(hh0[1]))
			chr.first = Gen::MATERNAL;
		
		if (Gen::is_haplotype<Gen::H0>(hh0[0]) && Gen::is_haplotype<Gen::H1>(hh0[1]))
			chr.first = Gen::PATERNAL;
		
		if (Gen::is_haplotype<Gen::H1>(hh1[0]) && Gen::is_haplotype<Gen::H0>(hh1[1]))
			chr.second = Gen::MATERNAL;
		
		if (Gen::is_haplotype<Gen::H0>(hh1[0]) && Gen::is_haplotype<Gen::H1>(hh1[1]))
			chr.second = Gen::PATERNAL;
		
		return chr;
	}
	
}


#endif /* IBD_hpp */

