//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Gen_hpp
#define Gen_hpp

#include <ctype.h>
#include <limits.h>
#include <stdint.h>

#include <algorithm>
#include <array>
#include <iomanip>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <set>
#include <sstream>


namespace Gen
{
	using value_t = uint8_t;
	using index_t = uint_fast8_t;
	
	using value_vector_t = std::vector<value_t>;
	using index_vector_t = std::vector<index_t>;
	
	
	static constexpr index_t ploidy = 2;
	static constexpr index_t hap_n = 3;
	static constexpr index_t gen_n = 4;
	
	
	using hap_t = value_t;
	using gen_t = value_t;
	
	using hap_vector_t = value_vector_t;
	using gen_vector_t = value_vector_t;
	
	using hap_count_t = std::array< size_t, hap_n >;
	using gen_count_t = std::array< size_t, gen_n >;
	
	using hap_pair_t = std::array<hap_t, ploidy>;
	using hap_pair_vector_t = std::vector<hap_t>;
	using hap_vector_pair_t = std::array<hap_pair_vector_t, ploidy>;
	
	
	typedef value_t compress_gen_t;
	typedef value_vector_t compress_gen_vector_t;
	
	
	// Haplotype index
	enum HapType : index_t
	{
		H0 = 0, // ref allele
		H1 = 1, // alt allele
		H_ = 2  // missing
	};
	
	// Genotype index
	enum GenType : index_t
	{
		G0 = 0,  // homozygous ref
		G1 = 1,  // heterozygous
		G2 = 2,  // homozygous alt
		G_ = 3   // missing
	};
	
	// Distinguish haplotypes
	enum ChrType : index_t
	{
		MATERNAL = 0,
		PATERNAL = 1,
		UNPHASED = 2,
		CHR_VOID = 3,
	};
	
	
	//
	// Haplotype functions
	//
	
	constexpr bool valid_as_haplotype(const char c)
	{
		return (c >= '0' && c <= '9');
	}
	
	
	constexpr hap_t make_haplotype(const char c)
	{
		return (valid_as_haplotype(c)) ? c - '0': 10;
	}
	
	
	constexpr char haplotype_to_char(const hap_t h)
	{
		return (h >= 0 && h <= 9) ? '0' + h: '.';
	}
	
	
	constexpr HapType haplotype_to_index(const hap_t h)
	{
		return
		(h == make_haplotype('0')) ? H0:
		(h == make_haplotype('1')) ? H1:
		H_;
	}
	
	
	constexpr hap_t index_to_haplotype(const HapType i)
	{
		return (i == H0) ? 0: (i == H1) ? 1: 10;
	}
	
	
	template < HapType H >
	constexpr bool is_haplotype(const hap_t h)
	{
		return (H == haplotype_to_index(h));
	}
	
	template < HapType H >
	constexpr bool is_haplotype(const HapType i)
	{
		return (H == i);
	}
	
	
	//
	// Genotype functions
	//
	
	constexpr bool valid_as_genotype(const char c0, const char c1)
	{
		return (valid_as_haplotype(c0) && valid_as_haplotype(c1));
	}
	
	
	constexpr gen_t make_genotype(const char c0, const char c1, const bool ph)
	{
		return ((11 * make_haplotype(c0)) + make_haplotype(c1)) + (121 * ph);
	}
	
	
	constexpr bool genotype_is_phased(const gen_t g)
	{
		return (g > 120);
	}
	
	
	constexpr gen_t unphase_genotype(const gen_t g)
	{
		return (genotype_is_phased(g)) ? g - 121: g;
	}
	
	
	constexpr hap_pair_t genotype_to_haplotypes(gen_t g)
	{
		return {{ hap_t((unphase_genotype(g) - (unphase_genotype(g) % 11)) / 11), hap_t(unphase_genotype(g) % 11) }};
	}
	
	
	constexpr hap_t genotype_to_haplotype0(gen_t g)
	{
		return hap_t((unphase_genotype(g) - (unphase_genotype(g) % 11)) / 11);
	}
	
	
	constexpr hap_t genotype_to_haplotype1(gen_t g)
	{
		return hap_t(unphase_genotype(g) % 11);
	}
	
	
	constexpr GenType genotype_to_index(gen_t g)
	{
		return
		(unphase_genotype(g) == make_genotype('0', '0', false)) ? G0:
		(unphase_genotype(g) == make_genotype('0', '1', false)) ? G1:
		(unphase_genotype(g) == make_genotype('1', '0', false)) ? G1:
		(unphase_genotype(g) == make_genotype('1', '1', false)) ? G2:
		G_;
	}
	
	
	constexpr gen_t index_to_genotype(const GenType gt)
	{
		return
		(gt == G0) ? make_genotype('0', '0', false):
		(gt == G1) ? make_genotype('0', '1', false):
		(gt == G2) ? make_genotype('1', '1', false):
		make_genotype('.', '.', false);
	}
	
	constexpr gen_t index_to_genotype(const index_t i)
	{
		return (i == G0) ? 0: (i == G1) ? 1: (i == G2) ? 12: 120;
	}
	
	
	template < GenType G >
	constexpr bool is_genotype(const gen_t g)
	{
		return (G == genotype_to_index(g));
	}
	
	template < GenType G >
	constexpr bool is_genotype(const GenType i)
	{
		return (G == i);
	}
	
	
	//
	// Genotype compression
	//
	
	enum CmprGenType : value_t
	{
		CG00 = 0,
		CG01 = 1,  CG01P = 2,
		CG0X = 3,  CG0XP = 4,
		CG10 = 5,  CG10P = 6,
		CG11 = 7,
		CG1X = 8,  CG1XP = 9,
		CGX0 = 10, CGX0P = 11,
		CGX1 = 12, CGX1P = 13,
		CGXX = 14,
		CG__ = 15 // undefined
	};
	
	
	static constexpr value_t make_compressed(const hap_t M, const hap_t P, const bool ph)
	{
		return
		(is_haplotype<H0>(M)) ? ((is_haplotype<H0>(P)) ? CG00: (is_haplotype<H1>(P)) ? ((ph) ? CG01P: CG01): ((ph) ? CG0XP: CG0X)):
		(is_haplotype<H1>(M)) ? ((is_haplotype<H0>(P)) ? ((ph) ? CG10P: CG10): (is_haplotype<H1>(P)) ? CG11: ((ph) ? CG1XP: CG1X)):
		(is_haplotype<H_>(M)) ? ((is_haplotype<H0>(P)) ? ((ph) ? CGX0P: CGX0): (is_haplotype<H1>(P)) ? ((ph) ? CGX1P: CGX1): CGXX):
		CG__;
	}
	
	static constexpr value_t make_compressed(const gen_t g)
	{
		return make_compressed(genotype_to_haplotype0(g), genotype_to_haplotype1(g), genotype_is_phased(g));
	}
	
	
	static constexpr gen_t unmake_compressed(const value_t cg)
	{
		return
		(CG00 == cg) ? make_genotype('0', '0', true):
		(CG01 == cg) ? make_genotype('0', '1', false): (CG01P == cg) ? make_genotype('0', '1', true):
		(CG0X == cg) ? make_genotype('0', '.', false): (CG0XP == cg) ? make_genotype('0', '.', true):
		(CG10 == cg) ? make_genotype('1', '0', false): (CG10P == cg) ? make_genotype('1', '0', true):
		(CG11 == cg) ? make_genotype('1', '1', true):
		(CG1X == cg) ? make_genotype('1', '.', false): (CG1XP == cg) ? make_genotype('1', '.', true):
		(CGX0 == cg) ? make_genotype('.', '0', false): (CGX0P == cg) ? make_genotype('.', '0', true):
		(CGX1 == cg) ? make_genotype('.', '1', false): (CGX1P == cg) ? make_genotype('.', '1', true):
		make_genotype('.', '.', true);
	}
	
	
	inline value_vector_t compress_genotype_vector(const gen_vector_t & v, const size_t size)
	{
		static constexpr value_t off = static_cast<value_t>(CHAR_BIT * sizeof(value_t) / 2);
		static constexpr value_t max = static_cast<value_t>((off << 2) - 1);
		
		value_vector_t out;
		value_t        num = 0;
		
		out.reserve(size);
		out.push_back(make_compressed(v.front()));
		
		for (size_t i = 1, k = 0; i < size; ++i)
		{
			const value_t c = make_compressed(v.at(i));
			
			if (num < max && out[k] == c)
			{
				++num;
				continue;
			}
			
			if (num > 0)
			{
				out[k] |= (num << off);
			}
			
			out.push_back(c);
			num = 0;
			++k;
		}
		
		if (num > 0)
		{
			out.back() |= (num << off);
		}
		
		return out;
	}
	
	
	inline gen_vector_t decompress_genotype_vector(const value_vector_t & v, const size_t size, const size_t full)
	{
		static constexpr value_t off = static_cast<value_t>(CHAR_BIT * sizeof(value_t) / 2);
		static constexpr value_t max = static_cast<value_t>((off << 2) - 1);
		static constexpr value_t inv = ~max;
		
		gen_vector_t   g(size);
		value_vector_t n(size);
		
		size_t length = size;
		
		for (size_t i = 0; i < size; ++i)
		{
			g[i] = unmake_compressed(v.at(i) & max);
			n[i] = (v.at(i) & inv) >> off;
			
			length += n[i];
		}
		
		if (length != full)
			return gen_vector_t();
		
		gen_vector_t out(length);
		unsigned     num = 0;
		
		for (size_t i = 0, k = 0; i < length; ++i)
		{
			out[i] = g[k];
			
			if (num++ == n[k])
			{
				num = 0;
				++k;
			}
		}
		
		return out;
	}
};


#endif /* Gen_hpp */

