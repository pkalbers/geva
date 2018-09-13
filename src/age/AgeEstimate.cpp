//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "AgeEstimate.hpp"


using namespace Age;



// construct

Estimate::Estimate(const Param::Data para, const Gen::Marker::Key & site)
: focal(site)
, param(para)
, shared(0)
, others(0)
, lower(decimal_nil)
, upper(decimal_nil)
//, share_done(false)
//, share_time(decimal_nil)
{
	this->logsum = decimal_vector_t(this->param->nt, decimal_nil);
}


// include density result

bool Estimate::include(CCF & ccf, const bool share, const Gamete::Pair & pair)
{
	const decimal_vector_t & prior = this->param->prior;
	const size_t n_times = this->param->nt;
	
	const decimal_t min = prior.at(0);
	const decimal_t max = prior.at(n_times - 1);
	
	if (!ccf.good || !ccf.pass)
	{
		return false;
	}
	
	if (ccf.q25 >= max || ccf.q50 <= min || ccf.q50 >= max || ccf.q75 <= min)
	{
		ccf.good = false;
		return false;
	}
	
	
	// calculate log sums over time vector
	
	for (size_t i = 0; i < n_times; ++i)
	{
		this->logsum[i] += std::log(ccf.d.at(i));
	}
	
	
	// add to estimate
	
	if (share)
	{
		++this->shared;
		this->lower += std::log(ccf.q50);
	}
	else
	{
		++this->others;
		this->upper += std::log(ccf.q50);
	}
	
	return true;
}


// return final age estimate

CLE Estimate::estimate()
{
	static constexpr decimal_t half = static_cast<decimal_t>(0.5);
	static constexpr decimal_t ci_lower = static_cast<decimal_t>(0.025);
	static constexpr decimal_t ci_upper = static_cast<decimal_t>(0.975);
	
	const decimal_vector_t & prior = this->param->prior;
	const size_t n_times = this->param->nt;
	
	
	CLE out;
	
	if (this->shared == 0 || this->others == 0)
	{
		return out;
	}
	
	out.n_shared = this->shared;
	out.n_others = this->others;
	
	
	// estimator
	
	out.lower = this->lower / this->shared; // = geometric mean
	out.upper = this->upper / this->others;
	
	out.estim = (out.lower + out.upper) / decimal_two;
	
	out.estim = std::exp(out.estim);
	out.lower = std::exp(out.lower);
	out.upper = std::exp(out.upper);
	
	
	// summary stats
	
	decimal_vector_t seq(n_times, decimal_nil);
	
	size_t arglogmax = n_times;
	decimal_t logmax = decimal_neg;
	decimal_t seqsum = decimal_nil;
	
	for (size_t i = 0; i < n_times; ++i)
	{
		if (logmax < this->logsum[i])
		{
			logmax = this->logsum[i];
			arglogmax = i;
		}
	}
	
	if (arglogmax == 0 || arglogmax >= n_times - 1)
	{
		return out;
	}
	
	for (size_t i = 0; i < n_times; ++i)
	{
		seq[i] = std::exp(this->logsum[i] - logmax);
		seqsum += seq[i];
	}
	
	for (size_t i = 0; i < n_times; ++i)
	{
		seq[i] /= seqsum;
	}
	
	
	decimal_vector_t cumsum(n_times, decimal_nil);
	
	decimal_t max = seq[0];
	decimal_t min = std::abs(seq[0] - half);
	
	size_t argmax = 0;
	size_t argmin = 0;
	
	cumsum[0] = seq[0];
	
	for (size_t i = 1; i < n_times; ++i)
	{
		cumsum[i] = cumsum[ i - 1 ] + seq[i];
		
		const decimal_t value = std::abs(cumsum[i] - half);
		
		if (max < seq[i])
		{
			max = seq[i];
			argmax = i;
		}
		
		if (min > value)
		{
			min = value;
			argmin = i;
		}
	}
	
	
	// mean
	
	out.mean = decimal_nil;
	
	for (size_t i = 0; i < n_times; ++i)
	{
		out.mean += prior.at(i) * seq[i];
	}
	
	
	// mode
	
	out.mode = prior.at(argmax);
	
	
	// median
	
	out.median = prior.at(argmin);
	
	
	// pseudo 95% CI
	
	size_t l0 = 0, l1 = 0;
	size_t u0 = 0, u1 = 0;
	
	size_t k = 0;
	
	for (; k < n_times; ++k)
	{
		if (cumsum[k] < ci_lower)
		{
			l0 = k;
		}
		else
		{
			l1 = k;
			break;
		}
	}
	
	for (; k < n_times; ++k)
	{
		if (cumsum[k] < ci_upper)
		{
			u0 = k;
		}
		else
		{
			u1 = k;
			break;
		}
	}
	
	out.ci_025 = approx<decimal_t>(ci_lower, cumsum.at(l0), cumsum.at(l1), prior.at(l0), prior.at(l1));
	out.ci_975 = approx<decimal_t>(ci_upper, cumsum.at(u0), cumsum.at(u1), prior.at(u0), prior.at(u1));
	
	
	// output sequence
	
	for (size_t i = 0; i < n_times; ++i)
	{
		if (seq[i] > decimal_err)
			seq[i] /= max;
	}
	
	out.grid = std::move(seq);
	
	out.good = true;
	
	return out;
}

