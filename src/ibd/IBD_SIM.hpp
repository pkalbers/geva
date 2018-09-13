//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef IBD_SIM_hpp
#define IBD_SIM_hpp

#include <stdio.h>

#include "Gen.hpp"
#include "GenMarker.hpp"
#include "GenSample.hpp"
#include "IBD.hpp"


namespace IBD
{
	// True IBD from simulations
	namespace SIM
	{
		static const std::string name = "Simulation results";
		static const std::string abbr = "SIM";
		
		struct Truth
		{
			Gen::Sample::Key::Pair pair;
			std::pair< Gen::ChrType, Gen::ChrType > chr;
			Segment segment;
			bool shared;
			
			using Vector = std::vector< Truth >;
			using Map = std::unordered_map< size_t, Vector >;
		};
		
		// Container
		class Result
		{
		public:
			
			using Data = std::shared_ptr< Result >;
			
			// include simulated result
			void set(const size_t &, const Truth &);
			
			// access map
			Truth::Map    const & get() const;
			Truth::Vector const & get(const Gen::Marker::Key &) const;
			
			// get number
			size_t size() const;
			
		private:
			
			Truth::Map map;
		};
	}
}

#endif /* IBD_SIM_hpp */
