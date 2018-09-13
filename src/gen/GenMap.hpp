//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef GenMap_hpp
#define GenMap_hpp

#include <cfloat>
#include <cmath>
#include <cstdlib>

#include "Approx.h"
#include "Decimal.h"

#include "Gen.hpp"


namespace Gen
{
	// genetic map container
	class Map
	{
	public:
		
		typedef int                               chromosome_t;
		typedef size_t                            position_t;
		typedef std::vector< chromosome_t >       chromosome_list_t;
		typedef chromosome_list_t::const_iterator chromosome_list_i;
		
		// return values
		struct Element
		{
			// construct
			Element(const decimal_t, const decimal_t);
			
			// check if valid
			bool valid() const;
			
			const decimal_t rate;
			const decimal_t dist;
		};
		
		
		// constructs
		Map();
		Map(const decimal_t &);
		
		// insert new mapped site
		void set(const chromosome_t, const position_t, const decimal_t, const decimal_t);
		
		// approximate at position in chromosome
		Element get(const chromosome_t, const position_t) const;
		
		// return number of positions (per chromosome)
		size_t size() const;
		size_t size(const chromosome_t) const;
		
		// return list of chromosomes
		chromosome_list_t chromosomes() const;
		
	private:
		
		typedef std::pair< chromosome_t, position_t> gmap_kt; // key type
		typedef Element                              gmap_mt; // mapped type
		typedef std::map< gmap_kt, gmap_mt >         gmap_t;  // map type
		typedef gmap_t::const_iterator               gmap_i;  // map iterator
		
		
		const decimal_t rec_rate; // constant recombination rate
		const bool      constant; // flag that rate is constant
		
		gmap_t map; // genetic map
		size_t sum; // number of all mapped sites
	};
}


#endif /* GenMap_hpp */

