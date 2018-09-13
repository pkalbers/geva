//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBD_FGT_hpp
#define IBD_FGT_hpp

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "IBD.hpp"

#include "Gen.hpp"
#include "GenMarker.hpp"
#include "GenVariant.hpp"


namespace IBD
{
	// Four-Gametes Test
	namespace FGT
	{
		static const std::string name = "Four gamete test (FGT)";
		static const std::string abbr = "FGT";
		
		
		// Gamete frequencies
		class Frequency
		{
		public:
			
			using Data = std::shared_ptr< Frequency >;
			
			
			// construct
			Frequency(const Gen::Marker::Key &, const Gen::Sample::Key::Vector &, const Gen::Grid::Data);
			
			// return frequencies
			size_t sample(const Gen::Marker::Key &) const;
			size_t subset(const Gen::Marker::Key &) const;
			
			
		private:
			
			using freq_vector_t = std::vector< size_t >;
			
			
			// identify focal haplotype
			static Gen::ChrType identify(const Gen::Variant::Vector::Data, const Gen::Marker::Key &);
			
			
			const Gen::Grid::Data grid;
			freq_vector_t         freq;
		};
		
		
		// Four-Gamates algorithm
		class Algorithm
		{
		public:
			
			static constexpr unsigned four_gam_n = 4;
			
			using four_gam_t = std::array< Gen::hap_t, four_gam_n >;
			
			
			// construct
			Algorithm(const Gen::Variant::Vector::Data, const Gen::Variant::Vector::Data); // classic
			Algorithm(const Gen::Variant::Vector::Data, const Gen::ChrType, const Gen::Variant::Vector::Data, const Gen::ChrType, const Frequency::Data); // frequency dependent
			
			// detect segment
			Segment detect(const Gen::Marker::Key &);
			
			// print to stream
			void print(const Gen::Sample::Key::Pair &, const Gen::Marker::Key &, std::ostream & = std::cout) const;
			static void print_header(std::ostream & = std::cout);
			
			
		private:
			
			using sharer_pair_t = std::pair< Gen::Variant::Vector::Data, Gen::Variant::Vector::Data >;
			using chromo_pair_t = std::pair< Gen::ChrType, Gen::ChrType >;
			
			// detection methods
			Segment detect_classic(const Gen::Marker::Key &);
			Segment detect_frqbase(const Gen::Marker::Key &);
			
			
			const sharer_pair_t   sharers;
			const chromo_pair_t   chromos;
			const Frequency::Data gametes;
			const size_t length;
		};
	}
}


#endif /* IBD_FGT_hpp */

