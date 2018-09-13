//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Age.hpp"


using namespace Gen;
using namespace IBD;
using namespace Age;




// construct

Gamete::Gamete()
: individual(std::numeric_limits<size_t>::max())
, chromosome(CHR_VOID)
{}

Gamete::Gamete(const Gen::Sample::Key & ind, const Gen::ChrType chr)
: individual(ind)
, chromosome(chr)
{}

Gamete::Gamete(const Gamete & other)
: individual(other.individual)
, chromosome(other.chromosome)
{}

Gamete::Gamete(Gamete && other)
: individual(other.individual)
, chromosome(other.chromosome)
{}


// assign

Gamete & Gamete::operator = (const Gamete & other)
{
	this->individual = other.individual;
	this->chromosome = other.chromosome;
	return *this;
}

Gamete & Gamete::operator = (Gamete && other)
{
	this->individual = other.individual;
	this->chromosome = other.chromosome;
	return *this;
}


// sort

bool Gamete::operator <  (const Gamete & other) const
{
	return
	(this->individual.value <  other.individual.value) ? true:
	(this->individual.value == other.individual.value) ? (this->chromosome < other.chromosome): false;
}

bool Gamete::operator  > (const Gamete & other) const
{
	return
	(this->individual.value  > other.individual.value) ? true:
	(this->individual.value == other.individual.value) ? (this->chromosome > other.chromosome): false;
}

bool Gamete::operator == (const Gamete & other) const
{
	return (this->individual.value == other.individual.value) && (this->chromosome == other.chromosome);
}

bool Gamete::operator != (const Gamete & other) const
{
	return (this->individual.value != other.individual.value) || (this->chromosome != other.chromosome);
}



// construct

CCF::CCF()
: shape(0)
, rate(0)
, q25(decimal_nil)
, q50(decimal_nil)
, q75(decimal_nil)
, good(false)
, pass(true)
{}

CCF::CCF(const CCF & other)
: shape(other.shape)
, rate(other.rate)
, q25(other.q25)
, q50(other.q25)
, q75(other.q25)
, d(other.d)
, good(other.good)
, pass(other.pass)
{}

CCF::CCF(CCF && other)
: shape(other.shape)
, rate(other.rate)
, q25(other.q25)
, q50(other.q25)
, q75(other.q25)
, d(std::move(other.d))
, good(other.good)
, pass(other.pass)
{}


// assign

CCF & CCF::operator = (const CCF & other)
{
	this->shape = other.shape;
	this->rate  = other.rate;
	this->q25   = other.q25;
	this->q50   = other.q50;
	this->q75   = other.q75;
	this->d     = other.d;
	this->good  = other.good;
	this->pass  = other.pass;
	return *this;
}

CCF & CCF::operator = (CCF && other)
{
	this->shape = other.shape;
	this->rate  = other.rate;
	this->q25   = other.q25;
	this->q50   = other.q50;
	this->q75   = other.q75;
	this->d     = std::move(other.d);
	this->good  = other.good;
	this->pass  = other.pass;
	return *this;
}



// construct

CLE::CLE()
: good(false)
{}



// construct

Param::Param(const Gen::Grid::Data grid, const size_t ne, const decimal_t  mr, const size_t n_times, const decimal_t max_time, const bool logged)
: nt(n_times)
, Ng(grid->sample_size()) // genotypes
, Nh(grid->sample_size() * 2) // haplotypes
, Nm(grid->marker_size()) // markers
, Ne(static_cast<decimal_t>(ne))
, Mr(mr)
, theta(static_cast<decimal_t>(4 * ne) * mr) // default theta
, theta_est(false)
, theta_man(false)
, prior(n_times, decimal_nil)
, log_prior(n_times, decimal_nil)
{
	const decimal_t FourNe100 = (4.0 * this->Ne) / decimal_100;
	
	
	// settings
	
	this->include_prior = true;
	//this->use_mut_clock = true;
	//this->use_rec_clock = true;
	this->use_post_prob = false;
	this->use_hard_brks = true;
	
	this->run_mut_clock = true;
	this->run_rec_clock = true;
	this->run_cmb_clock = true;
	
	//this->apply_filter_fixed   = false;
	//this->apply_filter_detect  = true;
	this->apply_nearest_neighb = true;
	this->relax_nearest_neighb = true;
	this->use_tree_consistency = true;
	
	this->limit_sharers = 100;
	this->outgroup_size = 100;
	
	this->limit_sharers_max = 1000; // equal to ~0.5 million concordant pairs
	this->outgroup_size_max = 1000; // equal to ~0.5 million discordant pairs
	
	this->breakpt_range = 1000;
	
	this->nearest_range = 5000;
	
	this->threads = 1;
	
	
	// variables
	
	const size_t n_sites = grid->marker_size();
	
	this->boundary[LHS] = grid->marker().front().index;
	this->boundary[RHS] = grid->marker().back().index;
	
	this->position.resize(n_sites, decimal_nil);
	this->distance.resize(n_sites, decimal_nil);
	this->frequency.resize(n_sites, decimal_nil);
	this->log_het.resize(n_sites, decimal_nil);
	this->log_hom.resize(n_sites, decimal_nil);
	this->cum_log_hom.resize(n_sites, decimal_nil);
	
	
	for (size_t i = 0; i < n_sites; ++i)
	{
		Marker const & marker = grid->marker(i);
		
		this->position[i] = static_cast<decimal_t>(marker.position);
		
		this->distance[i] = marker.gen_dist * FourNe100;
		
		this->frequency[i] = static_cast<decimal_t>(marker.hap_count[H1]) / static_cast<decimal_t>(this->Nh);
		
		this->log_het[i] = std::log(decimal_two * this->frequency[i] * (decimal_one - this->frequency[i]));
		
		this->log_hom[i] = std::log(std::pow(this->frequency[i], decimal_two) + std::pow(decimal_one - this->frequency[i], decimal_two));
	}
	
	this->cum_log_hom[0] = this->log_hom[0];
	
	for (size_t i = 1; i < n_sites; ++i)
	{
		this->cum_log_hom[i] = this->cum_log_hom[ i - 1 ] + this->log_hom[i];
	}
	
	
	// write coalescent times prior
	
	if (logged)
	{
		const decimal_t gen0 = std::log(decimal_err);
		const decimal_t step = (std::log(max_time) - gen0) / static_cast<decimal_t>(n_times - 1);
		
		this->prior[0] = decimal_err;
		
		for (size_t i = 1; i < n_times; ++i)
		{
			this->prior.at(i) = std::exp(std::log(this->prior.at(i - 1)) + step);
		}
	}
	else
	{
		const decimal_t step = max_time / static_cast<decimal_t>(n_times - 1);
		
		for (size_t i = 1; i < n_times; ++i)
		{
			this->prior.at(i) = this->prior.at(i - 1) + step;
		}
		
		this->prior[0] = step / decimal_two;
	}
	
	for (size_t i = 0; i < n_times; ++i)
	{
		this->log_prior.at(i) = std::log(this->prior.at(i));
	}
}


// set theta manually

void Param::set_theta(const double th)
{
	this->theta = static_cast<decimal_t>(th);
	this->theta_man = true;
	
	if (this->theta_est)
		throw std::invalid_argument("Theta cannot be set manually when also estimated");
}


// estimate theta (Watterson estimator)

void Param::est_theta(const Grid::Data grid)
{
	const decimal_t range = static_cast<decimal_t>(grid->marker().back().position - grid->marker().front().position);
	
	decimal_t alpha = decimal_nil;
	
	for (size_t i = 1; i < this->Nh; ++i)
	{
		alpha += decimal_one / static_cast<decimal_t>(i);
	}
	
	this->theta = static_cast<decimal_t>(grid->marker_size()) / (alpha * range);
	this->theta_est = true;
	
	if (this->theta_man)
		throw std::invalid_argument("Theta cannot be estimated when also set manually");
}

