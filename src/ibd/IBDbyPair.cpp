//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "IBDbyPair.hpp"


using namespace Gen;
using namespace IBD;
using namespace Pair;


//
// Target of detection
//

// constructs

Target::Target(const size_t & k, const Sample::Key::Pair & shared, const Marker::Key::Set & focals)
: fk(k)
, pair(shared)
, site(shuffle(focals))
, missing(decimal_one)
{}


// shuffle focal sites

Marker::Key::Vector Target::shuffle(const Marker::Key::Set & focals)
{
	Marker::Key::Vector out(focals.begin(), focals.end());
	std::shuffle(out.begin(), out.end(), random_generator);
	return out;
}


//
// Result container
//

// construct

Result::Result(const Share::Data share)
: length(0)
{
	Share::Index::Map::const_iterator map, end = share->get().cend();
	
	for (map = share->get().cbegin(); map != end; ++map)
	{
		const size_t fk = map->first;
		const Share::Index::Pairs & index = map->second.pairs;
		
		Share::Index::Pairs::const_iterator it, ti = index.cend();
		
		for (it = index.cbegin(); it != ti; ++it)
		{
			this->length += it->second.size();
			this->target.emplace_back( std::make_shared<Target>(fk, it->first, it->second) );
		}
	}
}


// return number of pairs

size_t Result::size() const
{
	return this->length;
}


// access list

Target::List const & Result::get() const
{
	return this->target;
}


// print to stream

size_t Result::print(std::ostream & stream) const
{
	size_t count = 0;
	
	// header
	stream << "MarkerID Fk SampleID0 SampleID1 Missing LHS RHS" << std::endl;
	
	// content
	
	Target::List::const_iterator it, ti = this->target.cend();
	
	for (it = this->target.cbegin(); it != ti; ++it)
	{
		const size_t fk = (*it)->fk;
		const Sample::Key::Pair & pair = (*it)->pair;
		const decimal_t miss = (*it)->missing;
		
		Target::Segmentation::const_iterator tract, tract_end = (*it)->map.cend();
		
		for (tract = (*it)->map.cbegin(); tract != tract_end; ++tract)
		{
			const Segment & segm = tract->first;
			
			Target::Segmentation::mapped_type::const_iterator site, site_end = tract->second.cend();
			
			for (site = tract->second.cbegin(); site != site_end; ++site)
			{
				stream << site->value << ' ';
				stream << fk << ' ';
				stream << pair.first.value << ' ';
				stream << pair.second.value << ' ';
				stream << std::fixed << std::setprecision(6) << miss << ' ';
				stream << segm[LHS] << ' ';
				stream << segm[RHS] << std::endl;
				
				++count;
			}
		}
	}
	
	return count;
}



//
// Exection of detection
//

// construct

Detect::Detect(const DetectMethod detectmethod, const decimal_t max_miss, const Target::Data target_pair, const Gen::Grid::Data grid, const HMM::Model::Data hmm_model)
: method(detectmethod)
, target(target_pair)
, source(grid)
, model(hmm_model)
, max_missing_rate(max_miss)
{
	if (this->method == DETECT_HMM && !this->model)
	{
		throw std::invalid_argument("HMM requires a model");
	}
}


// execute detection using selected methods

void Detect::run()
{
	const Variant::Vector::Data a = this->source->get(this->target->pair.first);
	const Variant::Vector::Data b = this->source->get(this->target->pair.second);
	
	this->target->missing = missing_rate(a->gen(), b->gen(), a->size());
	
	if (this->target->missing <= this->max_missing_rate)
	{
		switch (this->method)
		{
			case DETECT_DGT: this->dgt(a, b); break;
			case DETECT_FGT: this->fgt(a, b); break;
			case DETECT_HMM: this->hmm(a, b); break;
			case DETECT_SIM: throw std::logic_error("No idea why this happened!"); break;
			case DETECT_VOID: return;
		}
	}
}


// run methods

void Detect::dgt(const Variant::Vector::Data a, const Variant::Vector::Data b)
{
	DGT::Algorithm method(a, b);
	
	Marker::Key::Vector::const_iterator it, ti = this->target->site.cend();
	
	for (it = this->target->site.cbegin(); it != ti; ++it)
	{
		this->target->map[ method.detect(*it) ].insert(*it);
	}
}

void Detect::fgt(const Variant::Vector::Data a, const Variant::Vector::Data b)
{
	FGT::Algorithm method(a, b);
	
	Marker::Key::Vector::const_iterator it, ti = this->target->site.cend();
	
	for (it = this->target->site.cbegin(); it != ti; ++it)
	{
		this->target->map[ method.detect(*it) ].insert(*it);
	}
}

void Detect::hmm(const Variant::Vector::Data a, const Variant::Vector::Data b)
{
	throw std::string("Removed because redefined and not used anymore");
	//	Marker::Key::Vector::const_iterator it, ti = this->target->site.cend();
	//
	//	for (it = this->target->site.cbegin(); it != ti; ++it)
	//	{
	//		HMM::Algorithm method(a, b, this->model, false);
	//
	//		this->target->map[ method.detect(this->target->fk, *it) ].insert(*it);
	//	}
}

