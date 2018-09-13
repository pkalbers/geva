//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef AgeDensity_hpp
#define AgeDensity_hpp

#include <cmath>
#include <mutex>
#include <stdexcept>

#include "IBD.hpp"

#include "Age.hpp"


namespace Age
{
	// Cumulative density of coalescent times
	class Density
	{
	public:
		
		// construct
		Density(const Param::Data, const bool, const Gen::Marker::Key &, const IBD::Segment &);
		
		// include segment differences
		void differences(const IBD::SegDiff &);
		
		// include posterior probabilities
		void probability(const decimal_vector_t &, const decimal_vector_t &);
		void probability(decimal_vector_t &&, decimal_vector_t &&);
		
		// estimate cumulative coalescent function
		CCF estimate(const ClockType);
		
		
		const Gen::Marker::Key focal;
		const bool             share;
		const Param::Data      param;
		
	private:
		
		// execute main calculation
		CCF with_uncertainty(const ClockType); // include uncertainty through probability calculations
		CCF with_certainty(const ClockType);   // assume full certainty about detected breakpoints
		
		// calculate uncertainty
		void physical_distance();
		void genetic_distance();
		void approx_probability();
		
		// calculate log likelihood surface
		decimal_t likelihood_estimate(const ClockType, const IBD::LR, const decimal_t &) const;
		decimal_vector_t likelihood_surface(const ClockType);
		
		
		const IBD::Segment segment;
		IBD::SegDiff       segdiff;
		
		IBD::Side<size_t> length;
		
		IBD::Distribution gen_length;  // recombination (genetic) length
		IBD::Distribution gen_deltas;  // delta of rec. length
		IBD::Distribution phy_length;  // physical length
		IBD::Distribution prob_distr;  // breakpoint probabilities
		
		bool differences_flag; // segment differences provided
		bool probability_flag; // probabilities provided
	};
}


#endif /* AgeDensity_hpp */

