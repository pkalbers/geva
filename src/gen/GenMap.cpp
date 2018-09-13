//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "GenMap.hpp"


using namespace Gen;


// return values

// construct

Map::Element::Element(const decimal_t _rate, const decimal_t _dist)
: rate(_rate)
, dist(_dist)
{}


// check if valid

bool Map::Element::valid() const
{
	return !(this->rate < decimal_nil && this->dist < decimal_nil);
}



// genetic map container

// constructs

Map::Map()
: rec_rate(decimal_nil)
, constant(false)
{}

Map::Map(const decimal_t & rr)
: rec_rate(rr)
, constant(true)
{}


// insert new mapped site

void Map::set(const Map::chromosome_t chromosome, const Map::position_t position, const decimal_t recombrate, const decimal_t distance)
{
	if (this->constant)
	{
		throw std::runtime_error("Unable to set genetic map coordinates when a constant recombination rate was defined");
	}
	
	gmap_kt cp(chromosome, position);
	gmap_mt rd(recombrate, distance);
	
	// do not insert if already present
	if (this->map.count(cp) > 0)
	{
		throw std::string("Position duplicate for chromosome " + std::to_string(chromosome) + " at position " + std::to_string(position));
	}
	
	// insert
	std::pair<gmap_i, bool> ins = this->map.insert(std::make_pair(cp, rd));
	
	if (!ins.second)
	{
		throw std::runtime_error("Unexpected runtime error: Unable to read position into memory");
	}
	
	decimal_t lower = distance - 1.0;
	decimal_t upper = distance + 1.0;
	
	if (ins.first != this->map.cbegin())
	{
		lower = std::prev(ins.first)->second.dist;
	}
	
	if (std::next(ins.first) != this->map.cend())
	{
		upper = std::next(ins.first)->second.dist;
	}
	
	// ensure sorting of genomic distances
	if (lower <= distance && distance <= upper)
	{
		++this->sum;
		return;
	}
	
	// remove inserted
	this->map.erase(ins.first);
}


// approximate at position in chromosome

Map::Element Map::get(const Map::chromosome_t chromosome, const Map::position_t position) const
{
	static constexpr decimal_t centmega = static_cast<decimal_t>(1e8);
	
	// when constant rate
	if (this->constant)
	{
		return gmap_mt(this->rec_rate * centmega, this->rec_rate * static_cast<decimal_t>(position * 100));
	}
	
	
	// when variable rate as in genetic map file
	
	if (this->sum < 2)
	{
		return gmap_mt(0.0, 0.0);
	}
	
	gmap_kt cp(chromosome, position);
	
	gmap_i undef = this->map.cend();
	gmap_i first = this->map.cbegin();
	gmap_i point = this->map.lower_bound(cp);
	gmap_i lower, upper;
	
	bool point_is_chrom = (point != undef && point->first.first == chromosome);
	bool lower_is_chrom, upper_is_chrom;
	
	if (point_is_chrom)
	{
		// exact position
		if (point->first.second == position)
		{
			return point->second;
		}
		
		if (point != first)
		{
			upper = point;
			lower = std::prev(point);
			
			upper_is_chrom = point_is_chrom;
			lower_is_chrom = (lower != undef && lower->first.first == chromosome);
			
			// interpolate between two positions
			if (lower_is_chrom)
			{
				return gmap_mt(approx<decimal_t>(position, lower->first.second, upper->first.second, lower->second.rate, upper->second.rate),
							   approx<decimal_t>(position, lower->first.second, upper->first.second, lower->second.dist, upper->second.dist));
			}
		}
		lower = point;
		upper = std::next(point);
		
		lower_is_chrom = point_is_chrom;
		upper_is_chrom = (upper != undef && upper->first.first == chromosome);
		
		// before first position
		if (upper_is_chrom)
		{
			return gmap_mt(decimal_nil, lower->second.dist);
			
			//return gmap_mt(approx<decimal_t>(position, lower->first.second, upper->first.second, lower->second.rate, upper->second.rate),
			//			     approx<decimal_t>(position, lower->first.second, upper->first.second, lower->second.dist, upper->second.dist));
		}
	}
	else
	{
		if (point != first)
		{
			upper = std::prev(point);
			upper_is_chrom = (upper != undef && upper->first.first == chromosome);
			
			if (upper != first && upper_is_chrom)
			{
				lower = std::prev(upper);
				lower_is_chrom = (lower != undef && lower->first.first == chromosome);
				
				// after last position
				if (lower_is_chrom)
				{
					return gmap_mt(decimal_nil, upper->second.dist);
					
					//return gmap_mt(approx<decimal_t>(position, lower->first.second, upper->first.second, lower->second.rate, upper->second.rate),
					//			     approx<decimal_t>(position, lower->first.second, upper->first.second, lower->second.dist, upper->second.dist));
				}
			}
		}
	}
	
	return gmap_mt(0.0, 0.0);
}


// return number of positions (per chromosome)

size_t Map::size() const
{
	return this->sum;
}

size_t Map::size(const chromosome_t chromosome) const
{
	size_t n = 0;
	gmap_i i = this->map.lower_bound(gmap_kt(chromosome, 0));
	gmap_i j = this->map.cend();
	
	while (i != j && i->first.first == chromosome)
	{
		++n;
		++i;
	}
	
	return n;
}


// return list of chromosomes

Map::chromosome_list_t Map::chromosomes() const
{
	chromosome_list_t list;
	chromosome_t      last = 0;
	
	for (gmap_i i = this->map.cbegin(), z = this->map.cend(); i != z; ++i)
	{
		if (last != i->first.first)
		{
			list.push_back(i->first.first);
			last = i->first.first;
		}
	}
	
	return list;
}

