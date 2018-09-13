//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef GenMarker_hpp
#define GenMarker_hpp


#include <array>
#include <string>

#include "Decimal.h"
#include "Identity.h"

#include "Gen.hpp"


namespace Gen
{
	class Marker
	{
	private:
		
		// Alleles in marker
		class Allele
		{
		public:
			
			// constructs
			Allele();
			//Allele(const Allele &);
			//Allele(Allele &&);
			
			// parse allele string
			void parse(const std::string &);
			
			// return alleles
			std::string operator [] (const size_t) const;
			
			// check if bi-allelic SNP
			bool is_snp() const;
			
			// return number of alleles
			size_t size() const;
			
			// return raw string
			std::string str() const;
			
			
			static constexpr char sep = ','; // seperator used to distinguish alleles in string
			
		private:
			
			typedef int_least8_t             value_t;
			typedef std::vector<std::string> cache_t;
			
			
			// encode SNP to determine mask value
			static value_t mask_snp(char , char);
			
			value_t value; // snp mask index
			cache_t cache; // raw input string
		};
		
		
	public:
		
		using Key = Identity<Marker>;
		
		typedef std::vector< Marker >   Vector;
		typedef Vector::const_iterator  Iterator;
		typedef Vector::const_reference Reference;
		
		
		// constructs
		Marker();
		//Marker(const Marker &);
		//Marker(Marker &&);
		
		
		// count alleles and genotype
		void count(const gen_t);
		
		// print to stream
		void print(std::ostream & = std::cout) const;
		
		// return as string
		std::string str() const;
		
		
		Key         index;
		std::string label;
		
		int    chromosome;
		size_t position;
		
		Allele allele; // observed alleles
		
		hap_count_t hap_count;   // 0, 1, and else
		gen_count_t gen_count; // (0/0), (0/1 or 1/0), (1/1), and else
		
		decimal_t rec_rate; // recombination rate, from genetic map
		decimal_t gen_dist; // genetic distance, from genetic map
		
		
		static const std::string header; // header for printing
	};
}


#endif /* GenMarker_hpp */

