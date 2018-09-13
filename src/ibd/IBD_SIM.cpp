//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "IBD_SIM.hpp"


using namespace Gen;
using namespace IBD;
using namespace SIM;


// include simulated result

void Result::set(const size_t & site, const Truth & element)
{
	this->map[ site ].push_back(element);
}


// access map

Truth::Map const & Result::get() const
{
	return this->map;
}

Truth::Vector const & Result::get(const Gen::Marker::Key & key) const
{
	return this->map.at(key.value);
}

// get number
size_t Result::size() const
{
	size_t n = 0;
	
	Truth::Map::const_iterator it, ti = this->map.cend();
	
	for (it = this->map.cbegin(); it != ti; ++it)
	{
		n += it->second.size();
	}
	
	return n;
}
