//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "AgeDensity.hpp"


using namespace Gen;
using namespace IBD;
using namespace Age;



// construct

Density::Density(const Param::Data para, const bool shared, const Marker::Key & site, const IBD::Segment & segm)
: focal(site)
, share(shared)
, param(para)
, segment(segm)
, segdiff(0, 0)
, differences_flag(false)
, probability_flag(false)
{
	if (this->segment[LHS] > this->focal || this->focal > this->segment[RHS])
	{
		throw std::invalid_argument(std::string("Segment does not contain focal site") +
									std::string("\nLHS: ") + std::to_string(this->segment[LHS]) +
									std::string("\nRHS: ") + std::to_string(this->segment[RHS]) +
									std::string("\nFocal: ") + std::to_string(this->focal));
	}

	this->length[LHS] = std::min(this->param->breakpt_range, (this->focal.value - this->segment[LHS].value) + 1);
	this->length[RHS] = std::min(this->param->breakpt_range, (this->segment[RHS].value - this->focal.value) + 1);
}


// include segment differences

void Density::differences(const IBD::SegDiff & S)
{
	this->segdiff = S;

	this->differences_flag = true;
}


// include posterior probabilities

void Density::probability(const decimal_vector_t & _lhs, const decimal_vector_t & _rhs)
{
	decimal_vector_t lhs = _lhs;
	decimal_vector_t rhs = _rhs;

	this->probability(std::move(lhs), std::move(rhs));
}

void Density::probability(decimal_vector_t && lhs, decimal_vector_t && rhs)
{
	if (lhs.size() < this->length[LHS] || rhs.size() < this->length[RHS])
	{
		throw std::runtime_error("Provided posterior probabilities do not cover the whole segment");
	}


	// LHS

	this->prob_distr[LHS].resize(this->length[LHS], decimal_nil);

	for (size_t i = this->focal.value - this->segment[LHS].value, k = 0; k < this->length[LHS]; ++k, --i) // from focal site
	{
		this->prob_distr[LHS][k] = lhs.at(i); // logged already
	}


	// RHS

	this->prob_distr[RHS].resize(this->length[RHS], decimal_nil);

	for (size_t i = this->segment[RHS].value - this->focal.value, k = 0; k < this->length[RHS]; ++k, --i) // from focal site
	{
		this->prob_distr[RHS][k] = rhs.at(i); // logged already
	}


	this->probability_flag = true;
}


// estimate density

CCF Density::estimate(const ClockType clock)
{
	// exclude zero-size segments when focal allele is shared
//	if (this->share && this->length[LHS] == 1 && this->length[RHS] == 1)
//	{
//		return decimal_vector_t(0);
//	}

	if (this->param->use_hard_brks)
	{
		return this->with_certainty(clock);
	}
	return this->with_uncertainty(clock);
}


// execute main calculation

// include uncertainty through probability calculations

CCF Density::with_uncertainty(const ClockType clock)
{
	CCF ccf;

	ccf.d = likelihood_surface(clock);

	if (!this->share)
	{
		for (size_t i = 0, n = ccf.d.size(); i < n; ++i)
		{
			ccf.d[i] = decimal_one - ccf.d[i];
		}
	}

	return ccf;
}


// assume full certainty about detected breakpoints

// cumulative gamma function (Erlang)
inline decimal_t gamma_cdf(const size_t a, const decimal_t b, const decimal_t t)
{
	static constexpr decimal_t near_one = decimal_one - (decimal_err / 2);
	static constexpr decimal_t near_nil = (decimal_err / 2);

	const decimal_t tb = t * b;
	const decimal_t lg = std::log(tb);

	decimal_t cum = decimal_nil;

	for (size_t i = 0; i < a; ++i)
	{
		cum += std::exp(-1 * std::lgamma(i + decimal_one) + (lg * i) - tb);

		if (cum > near_one) return near_nil;
	}

	return (cum < near_nil) ? near_one: decimal_one - cum;

//	const decimal_t ntb  = -1.0 * t * b;
//	const decimal_t ltb = std::log(t * b);
//
//	decimal_t fac = decimal_nil;
//	decimal_t cum = decimal_nil;
//
//	for (size_t i = 0; i < a; ++i)
//	{
//		cum += std::exp(ntb + fac - std::lgamma(i + decimal_one));
//		fac += ltb;
//	}
//
//	return decimal_one - cum;
}

CCF Density::with_certainty(const ClockType clock)
{
	static constexpr decimal_t near_one = decimal_one - decimal_err;
	static constexpr decimal_t near_nil = decimal_err;
	static constexpr decimal_t q25 = static_cast<decimal_t>(0.25);
	static constexpr decimal_t q50 = static_cast<decimal_t>(0.50);
	static constexpr decimal_t q75 = static_cast<decimal_t>(0.75);


	const Segment & boundary = this->param->boundary;

	const decimal_vector_t & distance = this->param->distance;
	const decimal_vector_t & position = this->param->position;

	const decimal_vector_t & times = this->param->prior;
	const size_t  n_times = this->param->nt;
	const decimal_t theta = this->param->theta;

	//const bool incl_mut_clock = this->param->use_mut_clock;
	//const bool incl_rec_clock = this->param->use_rec_clock;


	// gamma cdf

	size_t    shape = 1;
	decimal_t rate  = decimal_one;

	if (clock == MUT_CLOCK || clock == CMB_CLOCK)
	{
		if (!this->share) // count focal site as difference when not shared
			shape += 1;

		decimal_t pos_lhs = 0;
		decimal_t pos_rhs = 0;
		
		if (this->segment[LHS].value != boundary[LHS])
			pos_lhs = static_cast<decimal_t>(position.at(this->segment[LHS].value) + position.at(this->segment[LHS].value + 1)) / decimal_two;
		else
			pos_lhs = static_cast<decimal_t>(position.at(this->segment[LHS].value) - 1);
		
		if (this->segment[RHS].value != boundary[RHS])
			pos_rhs = static_cast<decimal_t>(position.at(this->segment[RHS].value) + position.at(this->segment[RHS].value - 1)) / decimal_two;
		else
			pos_rhs = static_cast<decimal_t>(position.at(this->segment[RHS].value) + 1);

		shape += this->segdiff[LHS] + this->segdiff[RHS];
		rate  += std::abs(pos_rhs - pos_lhs) * theta;
	}

	if (clock == REC_CLOCK || clock == CMB_CLOCK)
	{
		decimal_t gen_lhs = 0;
		decimal_t gen_rhs = 0;
		
		if (this->segment[LHS].value != boundary[LHS])
		{
			shape += 1;
			gen_lhs = static_cast<decimal_t>(distance.at(this->segment[LHS].value) + distance.at(this->segment[LHS].value + 1)) / decimal_two;
		}
		else
			gen_lhs = static_cast<decimal_t>(distance.at(this->segment[LHS].value) - decimal_err);
			
		if (this->segment[RHS].value != boundary[RHS])
		{
			shape += 1;
			gen_rhs = static_cast<decimal_t>(distance.at(this->segment[RHS].value) + distance.at(this->segment[RHS].value + 1)) / decimal_two;
		}
		else
			gen_rhs = static_cast<decimal_t>(distance.at(this->segment[RHS].value) + decimal_err);

		rate  += std::abs(gen_rhs - gen_lhs) * decimal_two;
	}


	CCF ccf;

	ccf.shape = shape;
	ccf.rate  = rate;

	ccf.d = (this->share) ? decimal_vector_t(n_times, near_one): decimal_vector_t(n_times, near_nil);

	decimal_t m25 = decimal_one;
	decimal_t m50 = decimal_one;
	decimal_t m75 = decimal_one;
	size_t    i25 = 0;
	size_t    i50 = 0;
	size_t    i75 = 0;

	for (size_t i = 0; i < n_times; ++i)
	{
		const decimal_t value = gamma_cdf(shape, rate, times[i]);

		if (value > near_one)
			break;

		ccf.d[i] = (this->share) ? value: decimal_one - value;


		// Quantiles

		const decimal_t p25 = std::abs(value - q25);
		const decimal_t p50 = std::abs(value - q50);
		const decimal_t p75 = std::abs(value - q75);

		if (m25 > p25) { m25 = p25; i25 = i; }
		if (m50 > p50) { m50 = p50; i50 = i; }
		if (m75 > p75) { m75 = p75; i75 = i; }
	}

	ccf.q25 = times.at(i25);
	ccf.q50 = times.at(i50);
	ccf.q75 = times.at(i75);

	ccf.good = true;

	return ccf;
}


// calculate uncertainty

void Density::physical_distance()
{
	const Segment & boundary = this->param->boundary;

	const decimal_vector_t & position = this->param->position;

	const decimal_t focal_position = position.at(this->focal.value);

	this->phy_length[LHS].resize(this->length[LHS], decimal_nil);
	this->phy_length[RHS].resize(this->length[RHS], decimal_nil);


	// LHS
	{
		size_t site;

		const decimal_t break_position = (boundary[LHS] == this->segment[LHS]) ? position.at(this->segment[LHS].value): position.at(this->segment[LHS].value - 1);

		if (boundary[LHS] == this->segment[LHS])
		{
			this->phy_length[LHS][0] = (focal_position - break_position) + decimal_err;
		}
		else
		{
			this->phy_length[LHS][0] = (focal_position - ((position.at(this->segment[LHS]) + break_position) / decimal_two)) + decimal_err;
		}

		site = this->segment[LHS].value + 1;

		for (size_t k = 1; k < this->length[LHS]; ++k)
		{
			this->phy_length[LHS][k] = (focal_position - ((position.at(site) + position.at(site - 1)) / decimal_two)) + decimal_err;

			++site;
		}
	}


	// RHS
	{
		size_t site;

		const decimal_t break_position = (boundary[RHS] == this->segment[RHS]) ? position.at(this->segment[RHS].value): position.at(this->segment[RHS].value + 1);

		if (boundary[RHS] == this->segment[RHS])
		{
			this->phy_length[RHS][0] = (break_position - focal_position) + decimal_err;
		}
		else
		{
			this->phy_length[RHS][0] = (((break_position + position.at(this->segment[RHS])) / decimal_two) - focal_position) + decimal_err;
		}

		site = this->segment[RHS].value - 1;

		for (size_t k = 1; k < this->length[RHS]; ++k)
		{
			this->phy_length[RHS][k] = (((position.at(site + 1) + position.at(site)) / decimal_two) - focal_position) + decimal_err;

			--site;
		}
	}
}

void Density::genetic_distance()
{
	const Segment & boundary = this->param->boundary;

	const decimal_vector_t & distance = this->param->distance;

	const decimal_t focal_distance = distance.at(this->focal.value);

	this->gen_length[LHS].resize(this->length[LHS], decimal_nil);
	this->gen_length[RHS].resize(this->length[RHS], decimal_nil);

	this->gen_deltas[LHS].resize(this->length[LHS], decimal_nil);
	this->gen_deltas[RHS].resize(this->length[RHS], decimal_nil);


	// LHS
	{
		size_t site;

		const decimal_t break_distance = (boundary[LHS] == this->segment[LHS]) ? distance.at(this->segment[LHS].value): distance.at(this->segment[LHS].value - 1);

		site = this->segment[LHS].value;

		for (size_t k = 0; k < this->length[LHS]; ++k)
		{
			this->gen_length[LHS][k] = (focal_distance - distance.at(site)) + decimal_err;

			++site;
		}

		this->gen_deltas[LHS][0] = (distance.at(this->segment[LHS]) - break_distance) + decimal_err;

		site = this->segment[LHS].value + 1;

		for (size_t k = 1; k < this->length[LHS]; ++k)
		{
			this->gen_deltas[LHS][k] = (distance.at(site) - distance.at(site - 1)) + decimal_err;

			++site;
		}
	}


	// RHS
	{
		size_t site;

		const decimal_t break_distance = (boundary[RHS] == this->segment[RHS]) ? distance.at(this->segment[RHS].value): distance.at(this->segment[RHS].value + 1);

		site = this->segment[RHS].value;

		for (size_t k = 0; k < this->length[RHS]; ++k)
		{
			this->gen_length[RHS][k] = (distance.at(site) - focal_distance) + decimal_err;

			--site;
		}

		this->gen_deltas[RHS][0] = (break_distance - distance.at(this->segment[RHS])) + decimal_err;

		site = this->segment[RHS].value - 1;

		for (size_t k = 1; k < this->length[RHS]; ++k)
		{
			this->gen_deltas[RHS][k] = (distance.at(site + 1) - distance.at(site)) + decimal_err;

			--site;
		}
	}
}

void Density::approx_probability()
{
	if (this->probability_flag)
	{
		return;
	}

	const Segment & boundary = this->param->boundary;

	const decimal_vector_t & cum_log_hom = this->param->cum_log_hom;
	const decimal_vector_t & log_het     = this->param->log_het;

	this->prob_distr[LHS] = decimal_vector_t(this->length[LHS], decimal_nil);
	this->prob_distr[RHS] = decimal_vector_t(this->length[RHS], decimal_nil);


	const Side< decimal_t > brk_het((boundary[LHS] == this->segment[LHS]) ? decimal_nil: log_het.at(this->segment[LHS] - 1),
									(boundary[RHS] == this->segment[RHS]) ? decimal_nil: log_het.at(this->segment[RHS] + 1));


	// LHS
	{
		size_t site;

		if (boundary[LHS] == this->segment[LHS])
		{
			site = this->segment[LHS].value;

			for (size_t k = 1; k < this->length[LHS]; ++k)
			{
				this->prob_distr[LHS][k] = cum_log_hom.at(site);

				++site;
			}
		}
		else
		{
			site = this->segment[LHS].value - 1;

			const decimal_t brk_hom = cum_log_hom.at(site);

			for (size_t k = 0; k < this->length[LHS]; ++k)
			{
				this->prob_distr[LHS][k] = brk_het[LHS] + cum_log_hom.at(site) - brk_hom;

				++site;
			}
		}
	}


	// RHS
	{
		size_t site = this->segment[RHS].value;

		const decimal_t brk_hom = cum_log_hom.at(site);

		for (size_t k = 0; k < this->length[RHS]; ++k)
		{
			this->prob_distr[RHS][k] = brk_het[RHS] + brk_hom - cum_log_hom.at(site);

			--site;
		}
	}
}


// calculate log likelihood surface

// mutation clock
inline decimal_t mut_clock(const size_t & S, const decimal_t & L, const decimal_t & t, const decimal_t & theta)
{
	const decimal_t Ut = L * t * theta;

	return (S == 0) ? decimal_nil - Ut: (std::log(Ut) * S) - Ut;

	//return (S == 0) ? decimal_nil - (L * t * theta): (lt * S) - (L * t * theta);
}

// recombination clock
//inline decimal_t rec_clock(const decimal_t B, const decimal_t & D, const decimal_t & R, const decimal_t & t, const decimal_t & lt)
inline decimal_t rec_clock(const decimal_t & D, const decimal_t & R, const decimal_t & t)
{
	static constexpr decimal_t neg = static_cast<decimal_t>(-1);

	return std::log(decimal_one - std::exp(neg * D * t *0.5)) - (R * t *0.5);

	//return (B == 0) ? decimal_nil - (R * t): (std::log(t) * B) - (R * t);
}

// estimate likelihood at time
decimal_t Density::likelihood_estimate(const ClockType clock, const IBD::LR side, const decimal_t & time) const
{
	const decimal_t theta     = this->param->theta;
	//const Segment & boundary  = this->param->boundary;

	const bool incl_mut_clock = (clock == MUT_CLOCK || clock == CMB_CLOCK); //this->param->use_mut_clock;
	const bool incl_rec_clock = (clock == REC_CLOCK || clock == CMB_CLOCK); //this->param->use_rec_clock;

	const size_t n = this->length[side];
	const size_t S = this->segdiff[side];
	const decimal_vector_t & L = this->phy_length[side];
	const decimal_vector_t & D = this->gen_deltas[side];
	const decimal_vector_t & R = this->gen_length[side];
	const decimal_vector_t & P = this->prob_distr[side];

	//const decimal_t B = (boundary[side] == this->segment[side]) ? decimal_nil: decimal_one;

	decimal_vector_t d(n, decimal_nil);
	size_t argmax = n;
	decimal_t max(decimal_neg);
	decimal_t sum(decimal_nil);


	for (size_t k = 0; k < n; ++k)
	{
		if (P[k] < decimal_neg)
		{
			d[k] = P[k];
			continue;
		}

		const decimal_t mut = (incl_mut_clock) ? mut_clock(S, L[k], time, theta): decimal_nil;
		const decimal_t rec = (incl_rec_clock) ? rec_clock(D[k], R[k], time): decimal_nil;

		const decimal_t clk = mut + rec + P[k];

		//printf("%lu  %.8f  %.8f  %.8f\n", S, mut, rec, P[k]);

		if (max < clk)
		{
			max = clk;
			argmax = k;
		}

		d[k] = clk;
	}

	if (argmax == n)
	{
		return decimal_nil;
	}

	for (size_t k = 0; k < n; ++k)
	{
		sum += std::exp(d[k] - max);
	}

	return max + std::log(sum);
}


// calculate cumulative density
inline void cum_density(decimal_vector_t & d, const size_t n)
{
	size_t argmax = n;
	decimal_t max = decimal_neg;
	decimal_t sum = decimal_nil;

	for (size_t i = 0; i < n; ++i)
	{
		if (max < d[i])
		{
			max = d[i];
			argmax = i;
		}
	}

	if (argmax == n)
	{
		d = decimal_vector_t(0);
		return;
	}

	for (size_t i = 0; i < n; ++i)
	{
		d[i] = std::exp(d[i] - max);

		sum += d[i];
	}

	for (size_t i = 0; i < n; ++i)
	{
		d[i] /= sum;
	}

	for (size_t i = 1; i < n; ++i)
	{
		d[i] += d[i - 1];
	}


	// avoid hard 0 or 1
	for (size_t i = 0; i < n; ++i)
	{
		d[i] *= decimal_one - (decimal_two * decimal_err);
		d[i] += decimal_err;
	}
}

decimal_vector_t Density::likelihood_surface(const ClockType clock)
{
	const decimal_vector_t & time = this->param->prior;
	//const decimal_vector_t & logt = this->param->log_prior;
	const size_t  n_times = this->param->nt;
	const bool incl_prior = this->param->include_prior;

	if (clock == MUT_CLOCK || clock == CMB_CLOCK) //(this->param->use_mut_clock)
		this->physical_distance();

	if (clock == REC_CLOCK || clock == CMB_CLOCK) //(this->param->use_rec_clock)
		this->genetic_distance();

	this->approx_probability(); // unless provided


	decimal_vector_t llk(n_times, decimal_nil);


	for (size_t i = 0; i < n_times; ++i)
	{
		llk[i] += this->likelihood_estimate(clock, LHS, time[i]);
		llk[i] += this->likelihood_estimate(clock, RHS, time[i]);

		if (incl_prior)
		{
			llk[i] -= time[i];
		}
	}

	cum_density(llk, n_times);

	return llk;
}

