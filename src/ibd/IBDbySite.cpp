//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "IBDbySite.hpp"


using namespace IBD;
using namespace Site;
using namespace Gen;


//
// Target of detection
//

// constructs

Target::Target(const size_t & k, const Marker::Key & focal, const Sample::Key::Vector & shared)
: fk(k)
, focus(focal)
, share(shared)
{}



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
		const Share::Index::Sites & index = map->second.sites;
		
		Share::Index::Sites::const_iterator it, ti = index.cend();
		
		for (it = index.cbegin(); it != ti; ++it)
		{
			if (it->second.size() >= Share::minimum)
			{
				const size_t n = it->second.size();
				this->length += (n * (n - 1)) / 2;
				this->target.emplace_back( std::make_shared<Target>(fk, it->first, it->second) );
			}
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
		const Marker::Key & site = (*it)->focus;
		
		Target::Segmentation::const_iterator tract, tract_end = (*it)->map.cend();
		
		for (tract = (*it)->map.cbegin(); tract != tract_end; ++tract)
		{
			const Sample::Key::Pair & pair = tract->first;
			const Segment & segm = tract->second.first;
			const decimal_t miss = tract->second.second;
			
			stream << site.value << ' ';
			stream << fk << ' ';
			stream << pair.first.value << ' ';
			stream << pair.second.value << ' ';
			stream << std::fixed << std::setprecision(6) << miss << ' ';
			stream << segm[LHS] << ' ';
			stream << segm[RHS] << std::endl;
			
			++count;
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
	// walkabout each pair
	
	Sample::Key::Vector::const_iterator s0, s1, end0, end1;
	
	end1 = this->target->share.cend();
	end0 = std::prev(end1);
	
	for (s0 = this->target->share.cbegin(); s0 != end0; ++s0)
	{
		const Variant::Vector::Data a = this->source->get(*s0);
		
		for (s1 = std::next(s0); s1 != end1; ++s1)
		{
			const Variant::Vector::Data b = this->source->get(*s1);
			
			const Sample::Key::Pair pair(*s0, *s1);
			
			const decimal_t miss = missing_rate(a->gen(), b->gen(), a->size());
			
			if (miss < this->max_missing_rate)
			{
				switch (this->method)
				{
					case DETECT_DGT: this->dgt(pair, a, b, miss); break;
					case DETECT_FGT: this->fgt(pair, a, b, miss); break;
					case DETECT_HMM: this->hmm(pair, a, b, miss); break;
					case DETECT_SIM: throw std::logic_error("How did we get here?"); break;
					case DETECT_VOID: return;
				}
			}
		}
	}
}


// print to stream

void Detect::print(std::ostream & stream) const
{
	// walkabout each pair
	
	Sample::Key::Vector::const_iterator s0, s1, end0, end1;
	
	end1 = this->target->share.cend();
	end0 = std::prev(end1);
	
	for (s0 = this->target->share.cbegin(); s0 != end0; ++s0)
	{
		const Variant::Vector::Data a = this->source->get(*s0);
		
		for (s1 = std::next(s0); s1 != end1; ++s1)
		{
			const Variant::Vector::Data b = this->source->get(*s1);
			
			const Sample::Key::Pair pair(*s0, *s1);
			
			switch (this->method)
			{
				case DETECT_DGT:
				{
					DGT::Algorithm result(a, b);
					result.print(pair, this->target->focus, stream);
					break;
				}
				case DETECT_FGT:
				{
					FGT::Algorithm result(a, b);
					result.print(pair, this->target->focus, stream);
					break;
				}
				case DETECT_HMM:
				{
					throw std::string("Removed because redefined and not used anymore");
					//					HMM::Algorithm result(a, b, this->model, false);
					//					result.detect(this->target->fk, this->target->focus);
					//					result.posterior(false);
					//					result.print(pair, this->target->focus, stream);
					break;
				}
				case DETECT_SIM:
				{
					throw std::logic_error("Why is this happening?");
					break;
				}
				case DETECT_VOID: return;
			}
		}
	}
}


// run methods

void Detect::dgt(const Sample::Key::Pair & pair, const Variant::Vector::Data a, const Variant::Vector::Data b, const decimal_t & miss)
{
	DGT::Algorithm method(a, b);
	
	this->target->map[ pair ] = std::make_pair(method.detect(this->target->focus), miss);
}

void Detect::fgt(const Sample::Key::Pair & pair, const Variant::Vector::Data a, const Variant::Vector::Data b, const decimal_t & miss)
{
	FGT::Algorithm method(a, b);
	
	this->target->map[ pair ] = std::make_pair(method.detect(this->target->focus), miss);
}

void Detect::hmm(const Sample::Key::Pair & pair, const Variant::Vector::Data a, const Variant::Vector::Data b, const decimal_t & miss)
{
	throw std::string("Removed because redefined and not used anymore");
	//	HMM::Algorithm method(a, b, this->model, false);
	//
	//	this->target->map[ pair ] = std::make_pair(method.detect(this->target->fk, this->target->focus), miss);
}

