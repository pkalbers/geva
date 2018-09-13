//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBD_DGT_hpp
#define IBD_DGT_hpp

#include <iterator>
#include <set>
#include <stdexcept>
#include <string>

#include "IBD.hpp"

#include "GenMarker.hpp"
#include "GenVariant.hpp"


namespace IBD
{
	// Discordant homozygote genotype method
	namespace DGT
	{
		static const std::string name = "Discordant genotype test (DGT)";
		static const std::string abbr = "DGT";
		
		using discord_set_t = std::set< size_t >;
		
		
		// Allele frequencies
		class Frequency
		{
		public:
			
			using Data = std::shared_ptr< Frequency >;
			
			
			// construct
			Frequency(const Gen::Sample::Key::Vector &, const Gen::Grid::Data);
			
			// return frequencies
			size_t sample(const Gen::Marker::Key &) const;
			size_t subset(const Gen::Marker::Key &) const;
			
			
		private:
			
			using freq_vector_t = std::vector< size_t >;
			
			
			const Gen::Grid::Data grid;
			freq_vector_t         freq;
		};
		
		
		// Discordant homozygote genotype algorithm
		class Algorithm
		{
		public:
			
			// construct
			Algorithm(const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data); // classic
			Algorithm(const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data, const Frequency::Data); // frequency dependent
			
			
			// detect segment
			Segment detect(const Gen::Marker::Key &);
			
			// print to stream
			void print(const Gen::Sample::Key::Pair &, const Gen::Marker::Key &, std::ostream & = std::cout) const;
			static void print_header(std::ostream & = std::cout);
			
			
		private:
			
			//using discord_set_t = std::set< size_t >;
			using sharer_pair_t = std::pair< Gen::Variant::Vector::Data, Gen::Variant::Vector::Data >;
			
			
			// detection methods
			Segment detect_classic(const Gen::Marker::Key &);
			Segment detect_frqbase(const Gen::Marker::Key &);
			
			
			//discord_set_t brk;
			
			const sharer_pair_t   sharers;
			const Frequency::Data alleles;
			const size_t length;
		};
	}
}


#endif /* IBD_DGT_hpp */

