//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "IBD_HMM.hpp"


//#define DEBUG_HMM


using namespace Gen;
using namespace IBD;
using namespace HMM;


// Model probabilities

// construct

HMM::Model::Model(const size_t & effective_size, const size_t & sample_size, inits_list && inits_data_con, inits_list && inits_data_dis, emiss_list && emiss_data, dists_list && dists_data)
: Ne(effective_size)
, Nh(sample_size)
, do_iterative(false)
, inits_con(std::move(inits_data_con))
, inits_dis(std::move(inits_data_dis))
, emiss(std::move(emiss_data))
, dists(std::move(dists_data))
{}


// generate transitions for target

void HMM::Model::prepare_transition(const size_t & target)
{
	std::lock_guard<std::mutex> lock(this->guard);

	if (this->trans.count(target) != 0)
	{
		return;
	}

	const size_t size = this->dists.size();

	trans_list fk_trans(size);

	for (size_t i = 0; i < size; ++i)
	{
		fk_trans[i] = this->calc_trans_matrix(target, this->Ne, this->Nh, this->dists[i]);
	}

	this->trans[ target ] = std::move(fk_trans);
}


// return probabilities

HMM::Model::inits_list const & HMM::Model::initial_con()
{
	return this->inits_con;
}

HMM::Model::inits_list const & HMM::Model::initial_dis()
{
	return this->inits_dis;
}

HMM::Model::emiss_list const & HMM::Model::emission()
{
	return this->emiss;
}

HMM::Model::trans_list const & HMM::Model::transition(const size_t & target)
{
	std::lock_guard<std::mutex> lock(this->guard);

	return this->trans.at(target);
}


// return genetic distance

decimal_t HMM::Model::distance(const size_t & site) const
{
	return this->dists[site];
}


// calculate transition matrix

decimal_t HMM::Model::calc_expected_age(const size_t & fk, const size_t & nh)
{
	if (fk <=  1) return decimal_err;
	if (fk >= nh) return decimal_two;

	const decimal_t f = static_cast<decimal_t>(fk) / static_cast<decimal_t>(nh);

	return ((-1 * decimal_two * f) / (decimal_one - f)) * std::log(f);
}

HMM::Model::trans_type HMM::Model::calc_trans_matrix(const size_t & fk, const size_t & ne, const size_t & nh, const decimal_t & dist)
{
	static const decimal_t cent = static_cast<decimal_t>(100);

	const decimal_t xage = (fk == 0) ? decimal_one: calc_expected_age(fk, nh);
	const decimal_t yage = static_cast<decimal_t>(-4) * static_cast<decimal_t>(ne);
	const decimal_t prob = std::exp( xage * yage * dist / cent );
	//const decimal_t coal = 0.5 * (decimal_one - std::exp(-6.0 * xage));

	if (prob < decimal_nil)
		throw std::runtime_error("Transition probability is negative");

	if (prob > decimal_one)
		throw std::runtime_error("Transition probability is invalid");

	trans_type q;

	q[NON_STATE][NON_STATE] = decimal_one;
	q[NON_STATE][IBD_STATE] = decimal_nil;

	q[IBD_STATE][NON_STATE] = decimal_one - prob;
	q[IBD_STATE][IBD_STATE] = prob;

	//q[NON_STATE][IBD_STATE] = (coal / (decimal_one - coal)) * q[IBD_STATE][NON_STATE];
	//q[NON_STATE][NON_STATE] = decimal_one - q[NON_STATE][IBD_STATE];

	return q;
}



//
// Hidden Markov Model algorithm
//


// construct

Algorithm::Algorithm(const hap_vector_t a, const hap_vector_t b, const Model::Data m, const bool dis)
: obs(obs_hap_pair_vector(a, b))
, model(m)
, size(a.size())
, good(false)
, post(false)
, is_discord(dis)
{
	static constexpr size_t max = std::numeric_limits<size_t>::max();
	
	//	if (this->obs.size() == 0)
	//	{
	//		throw std::runtime_error("Observation sequence is empty");
	//	}
	
	//	if (random_coin())
	//	{
	//		this->focal_A = v0->hap(MATERNAL);
	//		this->other_A = v0->hap(PATERNAL);
	//	}
	//	else
	//	{
	//		this->focal_A = v0->hap(PATERNAL);
	//		this->other_A = v0->hap(MATERNAL);
	//	}
	//
	//	if (random_coin())
	//	{
	//		this->focal_B = v1->hap(MATERNAL);
	//		this->other_B = v1->hap(PATERNAL);
	//	}
	//	else
	//	{
	//		this->focal_B = v1->hap(PATERNAL);
	//		this->other_B = v1->hap(MATERNAL);
	//	}
	
	this->fk = max;
	this->focal = max;
}


// decode segment from Viterbi path

Segment Algorithm::detect(const size_t & fk_, const Marker::Key & site)
{
	distr_n length = size_distr(site, this->size);
	
	this->fk = fk_;
	this->focal = site;
	
	//
	//	// determine focal chromosomes
	//
	//	const hap_t sfA = this->focal_A.at(this->focal);
	//	const hap_t soA = this->other_A.at(this->focal);
	//	const hap_t sfB = this->focal_B.at(this->focal);
	//	const hap_t soB = this->other_B.at(this->focal);
	//
	//	if (is_haplotype<H0>(sfA) && is_haplotype<H1>(soA))
	//		this->focal_A.swap(this->other_A);
	//
	//	if (is_haplotype<H0>(sfB) && is_haplotype<H1>(soB))
	//		this->focal_B.swap(this->other_B);
	//
	//	this->obs = obs_hap_pair_vector(this->focal_A, this->focal_B);
	
	
	// check genotype at focal site
	
	const ObsHapPair focus = this->obs.at(this->focal.value);
	
	if (!is_obstype<H11>(focus) && !is_obstype<H01>(focus))
	{
		throw std::logic_error("Invalid haplotype at focal site");
	}
	
	
	// clear containers
	
	this->v_path = init_path_distr(length);
	this->v_prob = init_prob_distr(length);
	this->v_wght = init_distribution(length[LHS], length[RHS]);
	
	this->f_prob = init_prob_distr(length);
	this->f_wght = init_distribution(length[LHS], length[RHS]);
	
	this->b_prob = init_prob_distr(length);
	this->b_wght = init_distribution(length[LHS], length[RHS]);
	
	this->p_prob = init_prob_distr(length);
	
	
	// prepare model
	
	if (this->is_discord)
		this->model->prepare_transition(0);
	else
		this->model->prepare_transition(this->fk);
	
	
	// execute for each side
	
	this->execute_viterbi(length[LHS], LHS);
	this->execute_viterbi(length[RHS], RHS);
	
	
	// detect segment
	Side< size_t > dist(1, 1);
	
	// LHS
	for (size_t k = 1; k < length[LHS]; ++k)
	{
		dist[LHS] = k;
		if (this->v_path[LHS][k] == NON_STATE)
			break;
	}
	
	// RHS
	for (size_t k = 1; k < length[RHS]; ++k)
	{
		dist[RHS] = k;
		if (this->v_path[RHS][k] == NON_STATE)
			break;
	}
	
	return Segment(site.value - dist[LHS], site.value + dist[RHS]);
	
	
//	// switch
//	
//	if (this->v_path[LHS][0] == IBD_STATE)
//	{
//		this->execute_viterbi_switch(length[LHS], LHS);
//		
//		size_t x = 0;
//		
//		for (size_t k = 1; k < length[LHS]; ++k)
//		{
//			x = k;
//			if (this->v_path[LHS][k] == NON_STATE)
//				break;
//		}
//		
//		if (this->v_path[LHS][0] == IBD_STATE)
//		{
//			if (dist[LHS] > x)
//				dist[LHS] = x;
//		}
//	}
//	
//	if (this->v_path[RHS][0] == IBD_STATE)
//	{
//		this->execute_viterbi_switch(length[RHS], RHS);
//		
//		size_t x = 0;
//		
//		for (size_t k = 1; k < length[RHS]; ++k)
//		{
//			x = k;
//			if (this->v_path[RHS][k] == NON_STATE)
//				break;
//		}
//		
//		if (this->v_path[RHS][0] == IBD_STATE)
//		{
//			if (dist[RHS] > x)
//				dist[RHS] = x;
//		}
//	}
//	
//	return Segment(site.value - dist[LHS], site.value + dist[RHS]);
}


//Segment Algorithm::detect_iterative(const size_t & fk_, const Marker::Key & site)
//{
//	distr_n length = size_distr(site, this->size);
//
//	this->fk = fk_;
//	this->focal = site;
//
//
//	// determine focal chromosomes
//
//	const hap_t sfA = this->focal_A.at(this->focal);
//	const hap_t soA = this->other_A.at(this->focal);
//	const hap_t sfB = this->focal_B.at(this->focal);
//	const hap_t soB = this->other_B.at(this->focal);
//
//	if (is_haplotype<H0>(sfA) && is_haplotype<H1>(soA))
//		this->focal_A.swap(this->other_A);
//
//	if (is_haplotype<H0>(sfB) && is_haplotype<H1>(soB))
//		this->focal_B.swap(this->other_B);
//
//	this->obs = obs_hap_pair_vector(this->focal_A, this->focal_B);
//
//
//	// check genotype at focal site
//
//	const ObsHapPair focus = this->obs.at(this->focal.value);
//
//	if (!is_obstype<H11>(focus) && !is_obstype<H01>(focus))
//	{
//		throw std::logic_error("Invalid haplotype at focal site");
//	}
//
//
//	// clear containers
//
//	this->v_path = init_path_distr(length);
//	this->v_prob = init_prob_distr(length);
//	this->v_wght = init_distribution(length[LHS], length[RHS]);
//
//	this->f_prob = init_prob_distr(length);
//	this->f_wght = init_distribution(length[LHS], length[RHS]);
//
//	this->b_prob = init_prob_distr(length);
//	this->b_wght = init_distribution(length[LHS], length[RHS]);
//
//	this->p_prob = init_prob_distr(length);
//
//
//	// prepare model
//
//	if (this->is_discord)
//		this->model->prepare_transition(0);
//	else
//		this->model->prepare_transition(this->fk);
//
//
//	// execute for each side
//
//	this->execute_viterbi(length[LHS], LHS);
//	this->execute_viterbi(length[RHS], RHS);
//
//
//
//	hap_vector_t fA = this->focal_A;
//	hap_vector_t oA = this->other_A;
//	hap_vector_t fB = this->focal_B;
//	hap_vector_t oB = this->other_B;
//
//	// LHS
//
//	size_t L = 1;
//
//	for (; L < length[LHS]; ++L)
//	{
//		if (this->v_path[LHS][L] == NON_STATE)
//			break;
//	}
//
//	size_t lhs = L;
//
//	// iterate
//	try
//	{
//		while (lhs < length[LHS])
//		{
//			for (size_t k = L; k < length[LHS]; ++k)
//			{
//				const size_t i = site.value - k;
//				const hap_t  a = fA.at(i);
//				const hap_t  b = fB.at(i);
//
//				fA[i] = oA[i];
//				oA[i] = a;
//
//				fB[i] = oB[i];
//				oB[i] = b;
//			}
//
//			// A
//
//			size_t LA = L;
//
//			for (size_t k = L; k < length[LHS]; ++k)
//			{
//				const size_t i = site.value - k;
//				this->obs.at(i) = obs_hap_pair(fA[i], this->focal_B[i]);
//			}
//
//			this->execute_viterbi(length[LHS], LHS);
//
//			for (; LA < length[LHS]; ++LA)
//			{
//				if (this->v_path[LHS][LA] == NON_STATE)
//					break;
//			}
//
//			// B
//
//			size_t LB = L;
//
//			for (size_t k = L; k < length[LHS]; ++k)
//			{
//				const size_t i = site.value - k;
//				this->obs.at(i) = obs_hap_pair(this->focal_A[i], fB[i]);
//			}
//
//			this->execute_viterbi(length[LHS], LHS);
//
//			for (; LB < length[LHS]; ++LB)
//			{
//				if (this->v_path[LHS][LB] == NON_STATE)
//					break;
//			}
//
//			// A+B
//
//			size_t LC = L;
//
//			for (size_t k = L; k < length[LHS]; ++k)
//			{
//				const size_t i = site.value - k;
//				this->obs.at(i) = obs_hap_pair(fA[i], fB[i]);
//			}
//
//			this->execute_viterbi(length[LHS], LHS);
//
//			for (; LC < length[LHS]; ++LC)
//			{
//				if (this->v_path[LHS][LC] == NON_STATE)
//					break;
//			}
//
//
//
//			if (LA > L && LA > LB && LA > LC)
//			{
//				this->focal_A.swap(fA);
//				this->other_A.swap(oA);
//				lhs = LA;
//			}
//			else if (LB > L && LB > LA && LB > LC)
//			{
//				this->focal_B.swap(fB);
//				this->other_B.swap(oB);
//				lhs = LB;
//			}
//			else if (LC > L && LC > LB && LC > LA)
//			{
//				this->focal_A.swap(fA);
//				this->other_A.swap(oA);
//				this->focal_B.swap(fB);
//				this->other_B.swap(oB);
//				lhs = LC;
//			}
//			else
//			{
//				break;
//			}
//
//			for (size_t k = L; k < length[LHS]; ++k)
//			{
//				const size_t i = site.value - k;
//				this->obs.at(i) = obs_hap_pair(this->focal_A.at(i), this->focal_B.at(i));
//			}
//
//			L = lhs;
//		}
//	}
//	catch (const std::exception & ex)
//	{
//		throw std::logic_error(std::string("Iteration error [LHS]: ") + ex.what());
//	}
//
//
//
//	// RHS
//
//	size_t R = 1;
//
//	for (; R < length[RHS]; ++R)
//	{
//		if (this->v_path[RHS][R] == NON_STATE)
//			break;
//	}
//
//	size_t rhs = R;
//
//	// iterate
//	try
//	{
//		while (rhs < length[RHS])
//		{
//			for (size_t k = R; k < length[RHS]; ++k)
//			{
//				const size_t i = site.value + k;
//				const hap_t  a = fA.at(i);
//				const hap_t  b = fB.at(i);
//
//				fA[i] = oA[i];
//				oA[i] = a;
//
//				fB[i] = oB[i];
//				oB[i] = b;
//			}
//
//			// A
//
//			size_t RA = R;
//
//			for (size_t k = R; k < length[RHS]; ++k)
//			{
//				const size_t i = site.value + k;
//				this->obs.at(i) = obs_hap_pair(fA[i], this->focal_B[i]);
//			}
//
//			this->execute_viterbi(length[RHS], RHS);
//
//			for (; RA < length[RHS]; ++RA)
//			{
//				if (this->v_path[RHS][RA] == NON_STATE)
//					break;
//			}
//
//			// B
//
//			size_t RB = R;
//
//			for (size_t k = R; k < length[RHS]; ++k)
//			{
//				const size_t i = site.value + k;
//				this->obs.at(i) = obs_hap_pair(this->focal_A[i], fB[i]);
//			}
//
//			this->execute_viterbi(length[RHS], RHS);
//
//			for (; RB < length[RHS]; ++RB)
//			{
//				if (this->v_path[RHS][RB] == NON_STATE)
//					break;
//			}
//
//			// A+B
//
//			size_t RC = R;
//
//			for (size_t k = R; k < length[RHS]; ++k)
//			{
//				const size_t i = site.value + k;
//				this->obs.at(i) = obs_hap_pair(fA[i], fB[i]);
//			}
//
//			this->execute_viterbi(length[RHS], RHS);
//
//			for (; RC < length[RHS]; ++RC)
//			{
//				if (this->v_path[RHS][RC] == NON_STATE)
//					break;
//			}
//
//
//			if (RA > R && RA > RB && RA > RC)
//			{
//				this->focal_A.swap(fA);
//				this->other_A.swap(oA);
//				rhs = RA;
//			}
//			else if (RB > R && RB > RA && RB > RC)
//			{
//				this->focal_B.swap(fB);
//				this->other_B.swap(oB);
//				rhs = RB;
//			}
//			else if (RC > R && RC > RB && RC > RA)
//			{
//				this->focal_A.swap(fA);
//				this->other_A.swap(oA);
//				this->focal_B.swap(fB);
//				this->other_B.swap(oB);
//				rhs = RC;
//			}
//			else
//			{
//				break;
//			}
//
//			for (size_t k = R; k < length[RHS]; ++k)
//			{
//				const size_t i = site.value + k;
//				this->obs.at(i) = obs_hap_pair(this->focal_A.at(i), this->focal_B.at(i));
//			}
//
//			R = rhs;
//		}
//	}
//	catch (const std::exception & ex)
//	{
//		throw std::logic_error(std::string("Iteration error [RHS]: ") + ex.what());
//	}
//
//	return Segment((lhs == length[LHS]) ? site.value - (lhs - 1): site.value - lhs,
//				   (rhs == length[RHS]) ? site.value + (rhs - 1): site.value + rhs);
//}



// return Viterbi path decoding

prob_distr_t const & Algorithm::viterbi() const
{
	if (!this->good)
	{
		throw std::logic_error("Viterbi path probabilties were not calculated");
	}

	return this->v_prob;
}

Distribution const & Algorithm::weights() const
{
	if (!this->good)
	{
		throw std::logic_error("Viterbi path probabilties were not calculated");
	}

	return this->v_wght;
}

path_distr_t const & Algorithm::path() const
{
	if (!this->good)
	{
		throw std::logic_error("Viterbi path probabilties were not calculated");
	}

	return this->v_path;
}


// return posterior probabilities

prob_distr_t const & Algorithm::posterior(const bool logged)
{
	if (!this->good)
	{
		throw std::logic_error("Posterior probabilties were not calculated");
	}

	if (!this->post)
	{
		// execute for each side

		distr_n length = size_distr(this->focal, this->size);

		this->execute_posterior(length[LHS], LHS, logged);
		this->execute_posterior(length[RHS], RHS, logged);


		this->post = true;
	}

	return this->p_prob;
}

decimal_vector_t Algorithm::posterior(const LR side, const HiddenState state, const bool logged)
{
	prob_distr_t const & post_prob = this->posterior(logged);

	const size_t n = post_prob[side].size();

	decimal_vector_t p(n);

	for (size_t k = 0; k < n; ++k)
	{
		p[k] = post_prob[side][k][state];
	}

	return p;
}

Distribution const & Algorithm::fwd_weights() const
{
	if (!this->good)
	{
		throw std::logic_error("Posterior probabilties were not calculated");
	}

	return this->f_wght;
}

Distribution const & Algorithm::bwd_weights() const
{
	if (!this->good)
	{
		throw std::logic_error("Posterior probabilties were not calculated");
	}

	return this->b_wght;
}


// print to stream

void Algorithm::print(const Gen::Sample::Key::Pair & pair, const Gen::Marker::Key & site, std::ostream & stream) const
{
	Model::emiss_list const & Emiss = this->model->emission();

	distr_n length = size_distr(site, this->size);

	if (!this->good)
	{
		throw std::logic_error("Probabilties were not calculated");
	}

	if (!this->post)
	{
		throw std::logic_error("Posterior probabilties were not calculated");
	}


	// LHS

	for (size_t k = 0; k < length[LHS]; ++k)
	{
		const size_t i = site.value - k;
		const Model::emiss_type e = Emiss.at(i);
		const ObsHapPair o = this->obs[i];

		stream << pair.first << ' ' << pair.second << ' ' << site.value << ' ' << i << ' ';

		switch (this->v_path[LHS][k])
		{
			case NON_STATE:   stream << "NON "; break;
			case IBD_STATE:   stream << "IBD "; break;
			case STATE_UNDEF: stream << "--- "; break;
		}

		stream << std::fixed << std::setprecision(8) << e[NON_STATE][o] << ' ';
		stream << std::fixed << std::setprecision(8) << e[IBD_STATE][o] << ' ';

		const decimal_t pp_non = this->p_prob[LHS][k][NON_STATE];
		const decimal_t pp_ibd = this->p_prob[LHS][k][IBD_STATE];

		stream << std::fixed << std::setprecision(8) << ((pp_non < decimal_err) ? decimal_nil: (pp_non > decimal_one - decimal_err) ? decimal_one: pp_non) << ' ';
		stream << std::fixed << std::setprecision(8) << ((pp_ibd < decimal_err) ? decimal_nil: (pp_ibd > decimal_one - decimal_err) ? decimal_one: pp_ibd) << std::endl;
	}


	// RHS

	for (size_t k = 0; k < length[RHS]; ++k)
	{
		const size_t i = site.value + k;
		const Model::emiss_type e = Emiss.at(i);
		const ObsHapPair o = this->obs[i];

		stream << pair.first << ' ' << pair.second << ' ' << site.value << ' ' << i << ' ';

		switch (this->v_path[RHS][k])
		{
			case NON_STATE:   stream << "NON "; break;
			case IBD_STATE:   stream << "IBD "; break;
			case STATE_UNDEF: stream << "--- "; break;
		}

		stream << std::fixed << std::setprecision(8) << e[NON_STATE][o] << ' ';
		stream << std::fixed << std::setprecision(8) << e[IBD_STATE][o] << ' ';

		const decimal_t pp_non = this->p_prob[RHS][k][NON_STATE];
		const decimal_t pp_ibd = this->p_prob[RHS][k][IBD_STATE];

		stream << std::fixed << std::setprecision(8) << ((pp_non < decimal_err) ? decimal_nil: (pp_non > decimal_one - decimal_err) ? decimal_one: pp_non) << ' ';
		stream << std::fixed << std::setprecision(8) << ((pp_ibd < decimal_err) ? decimal_nil: (pp_ibd > decimal_one - decimal_err) ? decimal_one: pp_ibd) << std::endl;
	}
}

void Algorithm::print_header(std::ostream & stream)
{
	stream << "SampleID0 SampleID1 FocalID MarkerID ViterbiPath ObsProb_NON ObsProb_IBD PostProb_NON PostProb_IBD" << std::endl;
}



// Helper functions


inline void init_prob(prob_t & prob, const Model::inits_type & inits, const Model::emiss_type & emiss, const ObsHapPair obs)
{
#ifdef DEBUG_HMM
	if (is_obstype<H__>(obs)) // should never be called
	{
		//prob[NON_STATE] = inits[NON_STATE];
		//prob[IBD_STATE] = inits[IBD_STATE];
		throw std::logic_error("Unexpected haplotype at focal site");
	}
#endif

	prob[NON_STATE] = inits[NON_STATE] * emiss[NON_STATE][ obs ];
	prob[IBD_STATE] = inits[IBD_STATE] * emiss[IBD_STATE][ obs ];
	
//	prob[NON_STATE] = 0.5 * emiss[NON_STATE][ obs ]; // fixed prob
//	prob[IBD_STATE] = 0.5 * emiss[IBD_STATE][ obs ]; // fixed prob

#ifdef DEBUG_HMM
	if (prob[NON_STATE] < decimal_nil || prob[IBD_STATE] < decimal_nil)
	{
		throw std::runtime_error("Initial state probabilities equal zero");
	}
#endif
}


inline void init_prob_switch(prob_t & prob, const Model::inits_type & inits, const Model::emiss_type & emiss, const ObsHapPair obs)
{
#ifdef DEBUG_HMM
	if (is_obstype<H__>(obs)) // should never be called
	{
		//prob[NON_STATE] = inits[NON_STATE];
		//prob[IBD_STATE] = inits[IBD_STATE];
		throw std::logic_error("Unexpected haplotype at focal site");
	}
#endif
	
		prob[NON_STATE] = inits[NON_STATE] * emiss[IBD_STATE][ obs ]; // switched here
		prob[IBD_STATE] = inits[IBD_STATE] * emiss[NON_STATE][ obs ]; // switched here
	
//	prob[NON_STATE] = 0.5 * emiss[IBD_STATE][ obs ]; // switched here
//	prob[IBD_STATE] = 0.5 * emiss[NON_STATE][ obs ]; // switched here
	
#ifdef DEBUG_HMM
	if (prob[NON_STATE] < decimal_nil || prob[IBD_STATE] < decimal_nil)
	{
		throw std::runtime_error("Initial state probabilities equal zero");
	}
#endif
}


inline void path_max(prob_t & prob, const prob_t & prev, const Model::trans_type & trans, const Model::emiss_type & emiss, const ObsHapPair obs)
{
	if (is_obstype<H__>(obs))
	{
		prob[NON_STATE] = decimal_one;
		prob[IBD_STATE] = decimal_one;
	}
	else
	{
		prob[NON_STATE] = emiss[NON_STATE][ obs ];
		prob[IBD_STATE] = emiss[IBD_STATE][ obs ];
	}
#ifdef DEBUG_HMM
	if (prob[NON_STATE] < decimal_min || prob[IBD_STATE] < decimal_min)
	{
		throw std::runtime_error("Emission probabilities equal zero");
	}
#endif
	// NON

	const decimal_t trans_non_non = prev[NON_STATE] * trans[NON_STATE][NON_STATE];
	const decimal_t trans_ibd_non = prev[IBD_STATE] * trans[IBD_STATE][NON_STATE];

	prob[NON_STATE] *= (trans_non_non > trans_ibd_non) ? trans_non_non: trans_ibd_non;


	// IBD

	const decimal_t trans_non_ibd = prev[NON_STATE] * trans[NON_STATE][IBD_STATE];
	const decimal_t trans_ibd_ibd = prev[IBD_STATE] * trans[IBD_STATE][IBD_STATE];

	prob[IBD_STATE] *= (trans_non_ibd > trans_ibd_ibd) ? trans_non_ibd: trans_ibd_ibd;
}


inline void path_max_switch(prob_t & prob, const prob_t & prev, const Model::trans_type & trans, const Model::emiss_type & emiss, const ObsHapPair obs)
{
	if (is_obstype<H__>(obs))
	{
		prob[NON_STATE] = decimal_one;
		prob[IBD_STATE] = decimal_one;
	}
	else
	{
		prob[NON_STATE] = emiss[IBD_STATE][ obs ]; // switched here
		prob[IBD_STATE] = emiss[NON_STATE][ obs ]; // switched here
	}
#ifdef DEBUG_HMM
	if (prob[NON_STATE] < decimal_min || prob[IBD_STATE] < decimal_min)
	{
		throw std::runtime_error("Emission probabilities equal zero");
	}
#endif
	// NON
	
	const decimal_t trans_non_non = prev[NON_STATE] * trans[NON_STATE][NON_STATE];
	const decimal_t trans_ibd_non = prev[IBD_STATE] * trans[IBD_STATE][NON_STATE];
	
	prob[NON_STATE] *= (trans_non_non > trans_ibd_non) ? trans_non_non: trans_ibd_non;
	
	
	// IBD
	
	const decimal_t trans_non_ibd = prev[NON_STATE] * trans[NON_STATE][IBD_STATE];
	const decimal_t trans_ibd_ibd = prev[IBD_STATE] * trans[IBD_STATE][IBD_STATE];
	
	prob[IBD_STATE] *= (trans_non_ibd > trans_ibd_ibd) ? trans_non_ibd: trans_ibd_ibd;
}


inline HiddenState path_argmax(const prob_t & prob)
{
#ifdef DEBUG_HMM
	const decimal_t non = prob[NON_STATE];
	const decimal_t ibd = prob[IBD_STATE];

	if (non < decimal_min && ibd < decimal_min)
	{
		throw std::runtime_error("Inferred state probabilities equal zero");
	}

	return (non > ibd) ? NON_STATE: IBD_STATE;
#else
	return (prob[NON_STATE] > prob[IBD_STATE]) ? NON_STATE: IBD_STATE;
#endif
}


inline HiddenState path_argmax(const prob_t & prob, const Model::trans_type & trans, const HiddenState state)
{
#ifdef DEBUG_HMM
	const decimal_t non = prob[NON_STATE] * trans[NON_STATE][ state ];
	const decimal_t ibd = prob[IBD_STATE] * trans[IBD_STATE][ state ];

	if (non < decimal_min && ibd < decimal_min)
	{
		throw std::runtime_error("Inferred state probabilities equal zero in backtrack");
	}

	return (non > ibd) ? NON_STATE: IBD_STATE;
#else
	return ((prob[NON_STATE] * trans[NON_STATE][ state ]) > (prob[IBD_STATE] * trans[IBD_STATE][ state ])) ? NON_STATE: IBD_STATE;
#endif
}


inline void post_fwd(prob_t & prob, const prob_t & prev, const Model::trans_type & trans, const Model::emiss_type & emiss, const ObsHapPair obs)
{
	// NON
	prob[NON_STATE]  = prev[NON_STATE] * trans[NON_STATE][NON_STATE];
	prob[NON_STATE] += prev[IBD_STATE] * trans[IBD_STATE][NON_STATE];

	// IBD
	prob[IBD_STATE]  = prev[NON_STATE] * trans[NON_STATE][IBD_STATE];
	prob[IBD_STATE] += prev[IBD_STATE] * trans[IBD_STATE][IBD_STATE];

	if (!is_obstype<H__>(obs))
	{
		prob[NON_STATE] *= emiss[NON_STATE][ obs ];
		prob[IBD_STATE] *= emiss[IBD_STATE][ obs ];
	}
}


inline void post_bwd(prob_t & prob, const prob_t & next, const Model::trans_type & trans, const Model::emiss_type & emiss, const ObsHapPair obs)
{
	const decimal_t emiss_non = (is_obstype<H__>(obs)) ? next[NON_STATE]: next[NON_STATE] * emiss[NON_STATE][ obs ];
	const decimal_t emiss_ibd = (is_obstype<H__>(obs)) ? next[IBD_STATE]: next[IBD_STATE] * emiss[IBD_STATE][ obs ];

	// NON
	prob[NON_STATE]  = trans[NON_STATE][NON_STATE] * emiss_non;
	prob[NON_STATE] += trans[NON_STATE][IBD_STATE] * emiss_ibd;

	// IBD
	prob[IBD_STATE]  = trans[IBD_STATE][NON_STATE] * emiss_non;
	prob[IBD_STATE] += trans[IBD_STATE][IBD_STATE] * emiss_ibd;
}


inline decimal_t scale(prob_t & prob)
{
#ifdef DEBUG_HMM
	if (prob[NON_STATE] < decimal_min && prob[IBD_STATE] < decimal_min)
	{
		throw std::runtime_error("Invalid state probabilities");
	}
#endif
	decimal_t weight = decimal_nil;

	if (prob[NON_STATE] > prob[IBD_STATE])
	{
		weight = prob[NON_STATE];
		prob[NON_STATE]  = decimal_one;
		prob[IBD_STATE] /= weight;
	}
	else
	{
		weight = prob[IBD_STATE];
		prob[IBD_STATE]  = decimal_one;
		prob[NON_STATE] /= weight;
	}

	return weight;
}


// execute main algorithm

void Algorithm::execute_viterbi(const size_t & length, const LR side)
{
	const obs_vector_t & Obs = this->obs;


	// get model

	Model::inits_list const & Inits = (this->is_discord) ? this->model->initial_dis(): this->model->initial_con();
	Model::emiss_list const & Emiss = this->model->emission();
	Model::trans_list const & Trans = (this->is_discord) ? this->model->transition(0): this->model->transition(this->fk);


	// choose side

	path_vector_t &    path = (side == LHS) ? this->v_path[LHS]: this->v_path[RHS];
	prob_vector_t &    prob = (side == LHS) ? this->v_prob[LHS]: this->v_prob[RHS];
	decimal_vector_t & wght = (side == LHS) ? this->v_wght[LHS]: this->v_wght[RHS];


	size_t next = this->focal.value;
	size_t prev = this->focal.value;


	// forward

	init_prob(prob[0], Inits.at(this->focal.value), Emiss.at(this->focal.value), Obs.at(this->focal.value));

	wght[0] = scale(prob[0]);

	for (size_t k = 1; k < length; ++k)
	{
		switch (side)
		{
			case LHS:  --next;  path_max(prob[k], prob[ k - 1 ], Trans.at(next), Emiss.at(next), Obs.at(next));  break;
			case RHS:  ++next;  path_max(prob[k], prob[ k - 1 ], Trans.at(prev), Emiss.at(next), Obs.at(next));  break;
		}

		wght[k] = scale(prob[k]);

		prev = next;
	}


	// backward

	path[ length - 1 ] = path_argmax(prob[ length - 1 ]);

	for (size_t k = length - 1; k > 0; --k)
	{
		switch (side)
		{
			case LHS:  ++prev;  path[ k - 1 ] = path_argmax(prob[ k - 1 ], Trans.at(next), path[k]);  break;
			case RHS:  --prev;  path[ k - 1 ] = path_argmax(prob[ k - 1 ], Trans.at(prev), path[k]);  break;
		}

		next = prev;
	}
}


void Algorithm::execute_viterbi_switch(const size_t & length, const LR side)
{
	const obs_vector_t & Obs = this->obs;
	
	
	// get model
	
	Model::inits_list const & Inits = (this->is_discord) ? this->model->initial_dis(): this->model->initial_con();
	Model::emiss_list const & Emiss = this->model->emission();
	Model::trans_list const & Trans = (this->is_discord) ? this->model->transition(0): this->model->transition(this->fk);
	
	
	// choose side
	
	path_vector_t &    path = (side == LHS) ? this->v_path[LHS]: this->v_path[RHS];
	prob_vector_t &    prob = (side == LHS) ? this->v_prob[LHS]: this->v_prob[RHS];
	decimal_vector_t & wght = (side == LHS) ? this->v_wght[LHS]: this->v_wght[RHS];
	
	
	size_t next = this->focal.value;
	size_t prev = this->focal.value;
	
	
	// forward
	
	init_prob_switch(prob[0], Inits.at(this->focal.value), Emiss.at(this->focal.value), Obs.at(this->focal.value));
	
	wght[0] = scale(prob[0]);
	
	for (size_t k = 1; k < length; ++k)
	{
		switch (side)
		{
			case LHS:  --next;  path_max_switch(prob[k], prob[ k - 1 ], Trans.at(next), Emiss.at(next), Obs.at(next));  break;
			case RHS:  ++next;  path_max_switch(prob[k], prob[ k - 1 ], Trans.at(prev), Emiss.at(next), Obs.at(next));  break;
		}
		
		wght[k] = scale(prob[k]);
		
		prev = next;
	}
	
	
	// backward
	
	path[ length - 1 ] = path_argmax(prob[ length - 1 ]);
	
	for (size_t k = length - 1; k > 0; --k)
	{
		switch (side)
		{
			case LHS:  ++prev;  path[ k - 1 ] = path_argmax(prob[ k - 1 ], Trans.at(next), path[k]);  break;
			case RHS:  --prev;  path[ k - 1 ] = path_argmax(prob[ k - 1 ], Trans.at(prev), path[k]);  break;
		}
		
		next = prev;
	}
}


void Algorithm::execute_forward(const size_t & length, const LR side)
{
	const obs_vector_t & Obs = this->obs;


	// get model

	Model::inits_list const & Inits = (this->is_discord) ? this->model->initial_dis(): this->model->initial_con();
	Model::emiss_list const & Emiss = this->model->emission();
	Model::trans_list const & Trans = (this->is_discord) ? this->model->transition(0): this->model->transition(this->fk);


	// choose side

	prob_vector_t &    fwdp = (side == LHS) ? this->f_prob[LHS]: this->f_prob[RHS];
	decimal_vector_t & fwdw = (side == LHS) ? this->f_wght[LHS]: this->f_wght[RHS];


	size_t next = this->focal.value;
	size_t prev = this->focal.value;


	// forward

	init_prob(fwdp[0], Inits.at(this->focal.value), Emiss.at(focal.value), Obs.at(this->focal.value));

	fwdw[0] = scale(fwdp[0]);


	for (size_t k = 1; k < length; ++k)
	{
		switch (side)
		{
			case LHS:  --next;  post_fwd(fwdp[k], fwdp[ k - 1 ], Trans.at(next), Emiss.at(next), Obs.at(next));  break;
			case RHS:  ++next;  post_fwd(fwdp[k], fwdp[ k - 1 ], Trans.at(prev), Emiss.at(next), Obs.at(next));  break;
		}

		fwdw[k] = scale(fwdp[k]);

		prev = next;
	}
}


void Algorithm::execute_backward(const size_t & length, const LR side)
{
	const obs_vector_t & Obs = this->obs;


	// get model

	Model::emiss_list const & Emiss = this->model->emission();
	Model::trans_list const & Trans = this->model->transition(this->fk);


	// choose side

	prob_vector_t    & bwdp = (side == LHS) ? this->b_prob[LHS]: this->b_prob[RHS];
	decimal_vector_t & bwdw = (side == LHS) ? this->b_wght[LHS]: this->b_wght[RHS];


	size_t next = (side == LHS) ? 0: (this->focal.value + length) - 1;
	size_t prev = next;


	// backward

	bwdp[ length - 1 ][NON_STATE] = decimal_one;
	bwdp[ length - 1 ][IBD_STATE] = decimal_one;
	bwdw[ length - 1 ] = decimal_one;

	for (size_t k = length - 1; k > 0; --k)
	{
		switch (side)
		{
			case LHS:  ++prev;  post_bwd(bwdp[ k - 1 ], bwdp[k], Trans.at(next), Emiss.at(next), Obs.at(next));  break;
			case RHS:  --prev;  post_bwd(bwdp[ k - 1 ], bwdp[k], Trans.at(prev), Emiss.at(next), Obs.at(next));  break;
		}

		bwdw[ k - 1 ] = scale(bwdp[ k - 1 ]);

		next = prev;
	}
}


void Algorithm::execute_posterior(const size_t length, const LR side, const bool logged)
{
	const size_t l = length - 1;


	// execute fwd-bwd algorithm

	this->execute_forward(length, side);
	this->execute_backward(length, side);


	// choose side

	prob_vector_t & pp = this->p_prob[side];

	prob_vector_t & fp = this->f_prob[side];
	prob_vector_t & bp = this->b_prob[side];

	decimal_vector_t & fw = this->f_wght[side];
	decimal_vector_t & bw = this->b_wght[side];


	// cumulative sum of log weights

	fw[0] = std::log(fw.at(0));
	bw[l] = std::log(bw.at(l));

	for (size_t f = 1, b = l - 1; f < length; ++f, --b)
	{
		fw[f] = fw.at(f - 1) + std::log(fw.at(f));
		bw[b] = bw.at(b + 1) + std::log(bw.at(b));
	}


	// calculate posteriors

	const decimal_t non = fw[l] + std::log(fp[l][NON_STATE]);
	const decimal_t ibd = fw[l] + std::log(fp[l][IBD_STATE]);
	const decimal_t sum = non + std::log(decimal_one + std::exp(ibd - non));

	for (size_t k = 0; k < length; ++k)
	{
		pp[k][NON_STATE] = std::log(fp[k][NON_STATE] * bp[k][NON_STATE]) + fw[k] + bw[k] - sum;
		pp[k][IBD_STATE] = std::log(fp[k][IBD_STATE] * bp[k][IBD_STATE]) + fw[k] + bw[k] - sum;
	}

	if (!logged)
	{
		for (size_t k = 0; k < length; ++k)
		{
			pp[k][NON_STATE] = std::exp(pp[k][NON_STATE]);
			pp[k][IBD_STATE] = std::exp(pp[k][IBD_STATE]);
		}
	}


	this->post = true;
}
