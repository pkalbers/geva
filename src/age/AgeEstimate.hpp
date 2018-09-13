//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef AgeEstimate_hpp
#define AgeEstimate_hpp

#include "Approx.h"
#include "Decimal.h"

#include "Gen.hpp"
#include "GenMarker.hpp"

#include "Age.hpp"
#include "AgeDensity.hpp"


namespace Age
{
	// Point estimate
	class Estimate
	{
	public:
		
		// construct
		Estimate(const Param::Data, const Gen::Marker::Key &);
		
		// include density result
		bool include(CCF &, const bool, const Gamete::Pair &);
		
		// return final age estimate
		CLE estimate();
		
		
	private:
		
		const Gen::Marker::Key focal;
		const Param::Data      param;
		
		decimal_vector_t logsum;
		
		size_t shared;
		size_t others;
		
		decimal_t lower;
		decimal_t upper;
	};
}


#endif /* AgeEstimate_hpp */

