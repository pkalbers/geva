//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "AgeInfer.hpp"


using namespace Gen;
using namespace IBD;
using namespace Age;


// Pair of individuals/chromosomes

// construct

Pair::Pair(const Gamete::Pair _pair, const bool _sharing)
: pair(_pair)
, sharing(_sharing)
, missing(decimal_one)
, segment(0, 0)
, segdiff(-1, -1)
, done(false)
{}


// print to file

void Pair::print_header(std::ostream & stream, const Param::Data param, const bool full)
{
	stream << "MarkerID Clock";

	if (full)
	{
		decimal_vector_t::const_iterator time, time_end = param->prior.cend();

		for (time = param->prior.cbegin(); time != time_end; ++time)
		{
			stream << " t" << std::fixed << std::setprecision(8) << *time;
		}
	}
	else
	{
//		stream << " Fk SampleID0 Chr0 SampleID1 Chr1 Shared Missing Pass SegmentLHS SegmentRHS S_LHS S_RHS Shape Rate q25 q50 q75";
		stream << " SampleID0 Chr0 SampleID1 Chr1 Shared Pass SegmentLHS SegmentRHS Shape Rate";
	}

	stream << std::endl;
}

void Pair::print(std::ostream & stream, const Param::Data param, const bool full) const
{
	if (!this->done)
		throw std::runtime_error("Pairwise analysis not completed");

	// lock site pointer
	if (this->site.expired())
		throw std::runtime_error("Unexpected site pointer deletion");

	Site::Data site_ptr = this->site.lock();


	// print each clock

	for (int c = 0; c < n_clocks; ++c)
	{
		if (!this->ccf[c].good)
			continue;

		stream << site_ptr->focus.value << ' ';

		if (c == 0) stream << "M ";
		if (c == 1) stream << "R ";
		if (c == 2) stream << "J ";

		if (full)
		{
			decimal_vector_t::const_iterator est, end = this->ccf[c].d.cend();

			for (est = this->ccf[c].d.cbegin(); est != end; ++est)
			{
				stream << ' ' << std::fixed << std::setprecision(8) << *est;
			}
		}
		else
		{
			const ChrType chr0 = this->pair.first.chromosome;
			const ChrType chr1 = this->pair.second.chromosome;

//			stream << site_ptr->fk << ' ';

			stream << this->pair.first.individual.value << ' ';

			switch (chr0)
			{
				case MATERNAL: stream << '0' << ' '; break;
				case PATERNAL: stream << '1' << ' '; break;
				case UNPHASED: stream << '5' << ' '; break;
				case CHR_VOID: stream << '9' << ' '; break;
			}

			stream << this->pair.second.individual.value << ' ';

			switch (chr1)
			{
				case MATERNAL: stream << '0' << ' '; break;
				case PATERNAL: stream << '1' << ' '; break;
				case UNPHASED: stream << '5' << ' '; break;
				case CHR_VOID: stream << '9' << ' '; break;
			}

			stream << ((this->sharing) ? '1': '0') << ' ';

//			stream << std::fixed << std::setprecision(8) << this->missing << ' ';

			stream << ((this->ccf[c].pass) ? '1': '0') << ' ';

			stream << this->segment[LHS] << ' ';
			stream << this->segment[RHS] << ' ';

//			stream << this->segdiff[LHS] << ' ';
//			stream << this->segdiff[RHS] << ' ';

			stream << this->ccf[c].shape << ' ';
			stream << std::fixed << std::setprecision(8) << this->ccf[c].rate << ' ';
//			stream << std::fixed << std::setprecision(8) << this->ccf[c].q25 << ' ';
//			stream << std::fixed << std::setprecision(8) << this->ccf[c].q50 << ' ';
//			stream << std::fixed << std::setprecision(8) << this->ccf[c].q75 << ' ';
		}

		stream << std::endl;
	}
}




// Select nearest neighbours


// Rank

// construct

Near::Rank::Rank(const Gamete::Pair & pair_, const size_t dist_)
: pair(pair_)
, dist(dist_)
, rand(random_number())
{}

Near::Rank::Rank(const Near::Rank & other)
: pair(other.pair)
, dist(other.dist)
, rand(other.rand)
{}

Near::Rank::Rank(Near::Rank && other)
: pair(other.pair)
, dist(other.dist)
, rand(other.rand)
{}


// sort

bool Near::Rank::operator < (const Near::Rank & other) const
{
	return (this->dist < other.dist) || ( (this->dist == other.dist) && (this->rand < other.rand) );
}

bool Near::Rank::operator > (const Near::Rank & other) const
{
	return (this->dist > other.dist) || ( (this->dist == other.dist) && (this->rand > other.rand) );
}



// Chunk

// construct

Near::Chunk::Chunk(const Marker::Key & focus, const Gamete & chromo, const hap_vector_t & hap, const Param::Data param)
: chr(chromo)
, lhs(std::min(param->nearest_range, focus.value))
, rhs(std::min(param->nearest_range, param->Nm - focus.value - 1))
{
	const size_t nl = std::min(param->nearest_range, focus.value);
	const size_t nr = std::min(param->nearest_range, param->Nm - focus.value - 1);

	// LHS
	if (nl > 0)
	{
		for (size_t i = focus.value - 1, j = 0; j < nl; --i, ++j)
		{
			lhs.at(j) = hap.at(i);
		}
	}

	// RHS
	if (nr > 0)
	{
		for (size_t i = focus.value + 1, j = 0; j < nr; ++i, ++j)
		{
			rhs.at(j) = hap.at(i);
		}
	}
}


// get Hamming distance

size_t Near::Chunk::dist(const Near::Chunk & other) const ////////////////////////////////
{
	hap_vector_t::const_iterator l0, l1, l_end = this->lhs.cend();
	hap_vector_t::const_iterator r0, r1, r_end = this->rhs.cend();

	size_t d = 0;

	// LHS
	for (l0 = this->lhs.cbegin(), l1 = other.lhs.cbegin(); l0 != l_end; ++l0, ++l1)
	{
		if ((is_haplotype<H0>(*l0) && is_haplotype<H1>(*l1)) ||
			(is_haplotype<H1>(*l1) && is_haplotype<H0>(*l1)))
			++d;
	}

	// RHS
	for (r0 = this->rhs.cbegin(), r1 = other.rhs.cbegin(); r0 != r_end; ++r0, ++r1)
	{
		if ((is_haplotype<H0>(*r0) && is_haplotype<H1>(*r1)) ||
			(is_haplotype<H1>(*r1) && is_haplotype<H0>(*r1)))
			++d;
	}

	return d;
}


// Select nearest neighbours

// construct

Near::Near(const size_t & fk, const Marker::Key & focus, const Grid::Data grid, const Param::Data param)
: pool(param->threads, &Hold::run)
{
	this->ins.reserve(fk);
	this->out.reserve(param->Nh - fk);

	for (size_t i = 0; i < param->Ng; ++i)
	{
		this->pool.task(Hold(this, i, focus, grid, param));

		//		Chromo chr_mat;
		//		chr_mat.individual = i;
		//		chr_mat.chromosome = MATERNAL;
		//
		//		Chromo chr_pat;
		//		chr_pat.individual = i;
		//		chr_pat.chromosome = PATERNAL;
		//
		//
		//		const Variant::Vector::Data indv = grid->get(i);
		//
		//		const hap_vector_t hap_mat = indv->hap(MATERNAL);
		//		const hap_vector_t hap_pat = indv->hap(PATERNAL);
		//
		//
		//		if (is_haplotype<H1>(hap_mat.at(focus.value)))
		//		{
		//			this->ins.push_back( Chunk(focus, chr_mat, hap_mat, param) );
		//		}
		//		if (is_haplotype<H0>(hap_mat.at(focus.value)))
		//		{
		//			this->out.push_back( Chunk(focus, chr_mat, hap_mat, param) );
		//		}
		//
		//		if (is_haplotype<H1>(hap_pat.at(focus.value)))
		//		{
		//			this->ins.push_back( Chunk(focus, chr_pat, hap_pat, param) );
		//		}
		//		if (is_haplotype<H0>(hap_pat.at(focus.value)))
		//		{
		//			this->out.push_back( Chunk(focus, chr_pat, hap_pat, param) );
		//		}
	}

	this->pool.open();
	this->pool.exec();
	this->pool.wait();
}



// Near: Hold

Near::Hold::Hold(Near * _near, const Gen::Sample::Key & _key, const Marker::Key & focus, const Grid::Data grid, const Param::Data param)
: ptr(_near)
, key(_key)
, foc(focus)
, grd(grid)
, par(param)
{}

void Near::Hold::run()
{
	Gamete chr_mat(this->key, MATERNAL);
	Gamete chr_pat(this->key, PATERNAL);


	const Variant::Vector::Data indv = this->grd->get(this->key);

	const hap_vector_t hap_mat = indv->hap(MATERNAL);
	const hap_vector_t hap_pat = indv->hap(PATERNAL);

	Chunk chk_mat(this->foc, chr_mat, hap_mat, this->par);
	Chunk chk_pat(this->foc, chr_pat, hap_pat, this->par);


	std::lock_guard<std::mutex> guard(ptr->lock);


	if      (is_haplotype<H1>(hap_mat.at(this->foc.value)))
		ptr->ins.push_back( std::move(chk_mat) );
	else if (is_haplotype<H0>(hap_mat.at(this->foc.value)))
		ptr->out.push_back( std::move(chk_mat) );

	if      (is_haplotype<H1>(hap_pat.at(this->foc.value)))
		ptr->ins.push_back( std::move(chk_pat) );
	else if (is_haplotype<H0>(hap_pat.at(this->foc.value)))
		ptr->out.push_back( std::move(chk_pat) );
}


// perform all pairwise comparisons

bool Near::pairwise(const Param::Data param)
{
	size_t n_ins = this->ins.size();
	size_t n_out = this->out.size();

	if (n_ins < 2 || n_out < 2)
		return false;


	// concord

	Chunk::List::const_iterator c0, c0_end = std::prev(this->ins.cend());
	Chunk::List::const_iterator c1, c1_end = this->ins.cend();

	for (c0 = this->ins.cbegin(); c0 != c0_end; ++c0)
	{
		for (c1 = std::next(c0); c1 != c1_end; ++c1)
		{
			this->concord.push_back(Rank(std::make_pair(c0->chr, c1->chr), 0)); // to have random sorting
			//this->concord.push_back(Rank(std::make_pair(c0->chr, c1->chr), c0->dist(*c1)));
		}
	}


	// discord

	Chunk::List::const_iterator d0, d0_end = this->ins.cend();
	Chunk::List::const_iterator d1, d1_end = this->out.cend();

	for (d0 = this->ins.cbegin(); d0 != d0_end; ++d0)
	{
		for (d1 = this->out.cbegin(); d1 != d1_end; ++d1)
		{
			this->discord.push_back(Rank(std::make_pair(d0->chr, d1->chr), d0->dist(*d1)));
		}
	}


	// initial sorting

	this->concord.sort(); // randomly sorted
	this->discord.sort(); // prioritise lowest distance

	//this->concord.sort();
	//this->concord.reverse(); // prioritise highest distance


	// diversify sorting

	if (param->relax_nearest_neighb)
	{
		size_t num = 0;

		Rank::List dis;

		while (!this->discord.empty())
		{
			dis.push_back(this->discord.front());
			this->discord.pop_front();

			std::set<Gamete> unique;

			unique.insert(dis.back().pair.second);

			Rank::List::iterator it, end = this->discord.end();

			for (it = this->discord.begin(); it != end; )
			{
				if (unique.count(it->pair.second))
				{
					++it;
					continue;
				}

				unique.insert(it->pair.second);

				dis.push_back(*it);
				++num;
				it = this->discord.erase(it);
			}

			if (num > param->outgroup_size)
				break;
		}

		this->discord.swap(dis);


//		num = 0;
//
//		Rank::List con;
//
//		while (!this->concord.empty())
//		{
//			con.push_back(this->concord.front());
//			this->concord.pop_front();
//
//			std::set<Gamete> unique;
//
//			unique.insert(con.back().pair.first);
//			unique.insert(con.back().pair.second);
//
//			Rank::List::iterator it, end = this->concord.end();
//
//			for (it = this->concord.begin(); it != end; )
//			{
//				if (unique.count(it->pair.first) || unique.count(it->pair.second))
//				{
//					++it;
//					continue;
//				}
//
//				unique.insert(it->pair.first);
//				unique.insert(it->pair.second);
//
//				con.push_back(*it);
//				++num;
//				it = this->concord.erase(it);
//			}
//
//			if (num > param->limit_sharers)
//				break;
//		}
//
//		this->concord.swap(con);
	}

	return true;
}


// apply filter to pairs

// bool Near::filter(const Param::Data param)
// {
// 	if (!param->apply_filter_hamdist)
// 		return true;
//
//
// 	// filter concordant pairs
// 	{
// 		Chunk::List::const_iterator it, end = this->ins.cend();
//
// 		for (it = this->ins.cbegin(); it != end; ++it)
// 		{
// 			double c_avg = 0;
// 			size_t c_num = 0;
//
// 			double d_avg = 0;
// 			size_t d_num = 0;
//
// 			{
// 				Rank::List::const_iterator c, c_end = this->concord.cend();
// 				Rank::List::const_iterator d, d_end = this->discord.cend();
//
// 				for (c = this->concord.cbegin(); c != c_end; ++c)
// 				{
// 					if ((it->chr.individual == c->pair.first.individual  && it->chr.chromosome == c->pair.first.chromosome) ||
// 						(it->chr.individual == c->pair.second.individual && it->chr.chromosome == c->pair.second.chromosome))
// 					{
// 						c_avg += static_cast<double>(c->dist);
// 						++c_num;
// 					}
// 				}
//
// 				for (d = this->discord.cbegin(); d != d_end; ++d)
// 				{
// 					if (it->chr.individual == d->pair.first.individual && it->chr.chromosome == d->pair.first.chromosome)
// 					{
// 						d_avg += static_cast<double>(d->dist);
// 						++d_num;
// 					}
// 				}
// 			}
//
//
// 			c_avg /= static_cast<double>(c_num);
// 			d_avg /= static_cast<double>(d_num);
//
//
// 			if (c_avg > d_avg)
// 			{
// 				size_t size = 0;
//
// 				Rank::List::iterator c, c_end = this->concord.end();
// 				Rank::List::iterator d, d_end = this->discord.end();
//
// 				for (c = this->concord.begin(); c != c_end; )
// 				{
// 					if ((it->chr.individual == c->pair.first.individual  && it->chr.chromosome == c->pair.first.chromosome) ||
// 						(it->chr.individual == c->pair.second.individual && it->chr.chromosome == c->pair.second.chromosome))
// 					{
// 						c = this->concord.erase(c);
// 						continue;
// 					}
// 					++c;
//
// 					++size;
// 				}
//
// 				if (size == 0)
// 					return false;
//
// 				for (d = this->discord.begin(); d != d_end; )
// 				{
// 					if (it->chr.individual == d->pair.first.individual && it->chr.chromosome == d->pair.first.chromosome)
// 					{
// 						d = this->discord.erase(d);
// 						continue;
// 					}
// 					++d;
// 				}
// 			}
// 		}
// 	}
//
// 	return true;
// }



// SortNode

// construct

SortNode::SortNode(const decimal_t & _node, const Pair::List::iterator & _pair)
: node(_node)
, rand(random_number())
, pair(_pair)
{}

// sort

bool SortNode::operator < (const SortNode & other) const
{
	return (this->node < other.node) || ( (this->node == other.node) && (this->rand < other.rand) );
}

bool SortNode::operator > (const SortNode & other) const
{
	return (this->node > other.node) || ( (this->node == other.node) && (this->rand > other.rand) );
}



// Focal site shred by subset of individuals/chromosomes

// construct

Site::Site(const size_t & _fk, const Marker::Key & _focus, const Sample::Key::Vector & _share, const Grid::Data grid, const Param::Data param)
: fk(_fk)
, focus(_focus)
, share(_share)
, freq(nullptr)
, done(false)
{
	if (param->apply_nearest_neighb)
	{
		try
		{
			Near select(this->fk, this->focus, grid, param);

			if (select.pairwise(param))
			{
				size_t c_lim = 0;
				size_t d_lim = 0;

				Near::Rank::List::const_iterator c, c_end = select.concord.cend();
				Near::Rank::List::const_iterator d, d_end = select.discord.cend();

				for (c = select.concord.cbegin(); c != c_end; ++c)
				{
					this->list.push_back(std::make_shared< Pair >(c->pair, true));

					if (++c_lim == param->limit_sharers)
						break;
				}

				for (d = select.discord.cbegin(); d != d_end; ++d)
				{
					this->list.push_back(std::make_shared< Pair >(d->pair, false));

					if (++d_lim == param->outgroup_size)
						break;
				}

				if (c_lim == 0 || d_lim == 0)
					this->done = true;
			}
		}
		catch (const std::exception & ex)
		{
			throw std::runtime_error(std::string("Error while scanning for neaest neighbours at site: ") + std::to_string(this->focus.value) + std::string("\n") + + ex.what());
		}
	}
	else
	{
		const size_t n_sample = grid->sample_size();
		const size_t n_shared = this->share.size();
		const size_t n_others = n_sample - n_shared;
		const size_t n_outgrp = std::min(param->outgroup_size, n_others * n_shared);
		const std::unordered_set< size_t > sharer_set(this->share.cbegin(), this->share.cend());

		size_t n_paired = 0;


		// each sharer pair

		for (size_t i = 0; i < n_shared - 1; ++i)
		{
			for (size_t j = i + 1; j < n_shared; ++j)
			{
				Gamete::Pair pair;

				pair.first.individual = this->share[i];
				pair.first.chromosome = CHR_VOID; // to be specified during inference

				pair.second.individual = this->share[j];
				pair.second.chromosome = CHR_VOID; // to be specified during inference

				this->list.push_back(std::make_shared< Pair >(pair, true));

				++n_paired;
			}
		}


		// limit sharer pairs

		if (n_paired > param->limit_sharers)
		{
			Pair::List sub_list;

			random_vector_t sub_select = random_vector(param->limit_sharers, n_paired, true, true);

			random_vector_t::const_iterator it, ti = sub_select.cend();

			for (it = sub_select.cbegin(); it != ti; ++it)
			{
				sub_list.push_back(std::move(this->list[ *it ]));
			}

			this->list.swap(sub_list);

			n_paired = param->limit_sharers;
		}



		// subsample of outgroup pairs

		std::set< Sample::Key::Pair > unique;

		while (unique.size() < n_outgrp)
		{
			const size_t sharer = this->share.at(random_number(n_shared)).value;

			size_t other = random_number(n_sample);

			while (sharer_set.count(other) != 0)
			{
				other = random_number(n_sample);
			}

			unique.emplace(sharer, other);
		}

		std::set< Sample::Key::Pair >::iterator it, ti = unique.end();

		for (it = unique.begin(); it != ti; ++it)
		{
			Gamete::Pair pair;

			pair.first.individual = it->first;
			pair.first.chromosome = CHR_VOID; // to be specified during inference

			pair.second.individual = it->second;
			pair.second.chromosome = CHR_VOID; // to be specified during inference

			this->list.push_back(std::make_shared< Pair >(pair, false));
		}
	}
}


// construct for simulated results

Site::Site(const Gen::Marker::Key & site, const IBD::SIM::Result::Data simres)
: fk(0)
, focus(site)
, share(0)
, freq(nullptr)
, done(false)
{
	IBD::SIM::Truth::Vector const & index = simres->get(site.value);
	IBD::SIM::Truth::Vector::const_iterator it, ti = index.cend();

	for (it = index.cbegin(); it != ti; ++it)
	{
		Gamete::Pair pair;

		pair.first.individual = it->pair.first;
		pair.first.chromosome = it->chr.first;

		pair.second.individual = it->pair.second;
		pair.second.chromosome = it->chr.second;

		Pair::Data result = std::make_shared< Pair >(pair, it->shared);

		result->segment = it->segment;

		this->list.push_back(result);
	}
}


// calculate allele frequencies

void Site::frequency(const Grid::Data grid)
{
	std::lock_guard<std::mutex> lock(this->guard);

	if (!this->freq)
	{
		this->freq = std::make_shared< FGT::Frequency >(this->focus, this->share, grid);
	}
}


// estimate age

void Site::estimate(const Param::Data param)
{
	std::lock_guard<std::mutex> lock(this->guard);

	if (this->done)
	{
		return;
	}

	Pair::List::const_iterator it, ti = this->list.cend();

	for (it = this->list.cbegin(); it != ti; ++it)
	{
		if (!(*it)->done)
		{
			return;
		}
	}

	this->done = true;


	// raw estimate

	for (int c = 0; c < n_clocks; ++c)
	{
		const ClockType clock = static_cast<ClockType>(c);

		Estimate age(param, this->focus);

		size_t count = 0;

		for (it = this->list.cbegin(); it != ti; ++it)
		{
			if (age.include((*it)->ccf[clock], (*it)->sharing, (*it)->pair))
			{
				++count;
			}
		}

		if (count > 0)
		{
			this->raw[clock] = age.estimate();
		}
	}


	// adjusted estimate

	for (int c = 0; c < n_clocks; ++c)
	{
		const ClockType clock = static_cast<ClockType>(c);

		//
		this->filter(clock, param); // filter pairs
		//

		Estimate age(param, this->focus);

		size_t count = 0;

		for (it = this->list.cbegin(); it != ti; ++it)
		{
			if (age.include((*it)->ccf[clock], (*it)->sharing, (*it)->pair))
			{
				++count;
			}
		}

		if (count > 0)
		{
			this->adj[clock] = age.estimate();
		}
	}
}


// filter pairs based on CCF summary metric

void Site::filter(const ClockType clock, const Param::Data param)
{
	static constexpr decimal_t max_exclude_prop = static_cast<decimal_t>(0.5);

	const size_t n_times = param->nt;
	size_t ncon = 0;
	size_t ndis = 0;

	MinExclude::List minx(n_times);

	Pair::List::iterator it, ti = this->list.end();


	// count pairs

	for (it = this->list.begin(); it != ti; ++it)
	{
		if ((*it)->ccf[clock].good)
		{
			if ((*it)->sharing)
				++ncon;
			else
				++ndis;
		}
	}


	// fill exclusion list

	for (size_t i = 0; i < n_times; ++i)
	{
		MinExclude & x = minx[i];

		x.time = param->prior.at(i);
		x.ncon = 0;
		x.ndis = 0;
		x.wsum = decimal_nil;

		for (it = this->list.begin(); it != ti; ++it)
		{
			if (!(*it)->ccf[clock].good)
				continue;

			if ((*it)->sharing)
			{
				if ((*it)->ccf[clock].q50 > x.time) // median of CCF
					++x.ncon;
			}
			else
			{
				if ((*it)->ccf[clock].q50 < x.time) // median of CCF
					++x.ndis;
			}
		}

		x.wsum += static_cast<decimal_t>(x.ncon) / static_cast<decimal_t>(ncon);
		x.wsum += static_cast<decimal_t>(x.ndis) / static_cast<decimal_t>(ndis);
	}


	// determine threshold

	decimal_t min = decimal_max;
	size_t    arg = n_times;

	for (size_t i = 0; i < n_times; ++i)
	{
		if (min > minx[i].wsum)
		{
			min = minx[i].wsum;
			arg = i;
		}
	}

	if (arg == n_times)
		return;


	// find pairs below/above threshold

	const decimal_t cutoff = minx[arg].time;

	SortNode::List tcon, tdis;

	for (it = this->list.begin(); it != ti; ++it)
	{
		if (!(*it)->ccf[clock].good)
			continue;

		if ((*it)->sharing)
		{
			if ((*it)->ccf[clock].q50 > cutoff) // median of CCF
				tcon.push_back( SortNode((*it)->ccf[clock].q50, it) );
		}
		else
		{
			if ((*it)->ccf[clock].q50 < cutoff) // median of CCF
				tdis.push_back( SortNode((*it)->ccf[clock].q50, it) );
		}
	}


	// apply filter

	const size_t ntcon = tcon.size();
	const size_t ntdis = tdis.size();

	if (ncon > 1 && ntcon > 0)
	{
		const size_t max = std::floor(static_cast<decimal_t>(ncon) * max_exclude_prop);

		if (max > ntcon)
		{
			for (SortNode::List::const_reference ex : tcon)
			{
				(*ex.pair)->ccf[clock].pass = false;
			}
		}
		else
		{
			tcon.sort();
			tcon.reverse(); // sort highest first

			size_t num = 0;

			for (SortNode::List::const_reference ex : tcon)
			{
				(*ex.pair)->ccf[clock].pass = false;

				if (++num == max)
					break;
			}
		}
	}

	if (ndis > 1 && ntdis > 0)
	{
		const size_t max = std::floor(static_cast<decimal_t>(ndis) * max_exclude_prop);

		if (max > ntdis)
		{
			for (SortNode::List::const_reference ex : tdis)
			{
				(*ex.pair)->ccf[clock].pass = false;
			}
		}
		else
		{
			tdis.sort(); // sort lowest first

			size_t num = 0;

			for (SortNode::List::const_reference ex : tdis)
			{
				(*ex.pair)->ccf[clock].pass = false;

				if (++num == max)
					break;
			}
		}
	}
}


//void Site::filter_fixed(const ClockType clock, const Param::Data param)
//{
//	static constexpr decimal_t exclude_prop = static_cast<decimal_t>(0.1);
//
//	size_t ncon = 0;
//	size_t ndis = 0;
//
//	Pair::List::iterator it, ti = this->list.end();
//
//
//	// count pairs
//
//	for (it = this->list.begin(); it != ti; ++it)
//	{
//		if ((*it)->ccf[clock].good)
//		{
//			if ((*it)->sharing)
//				++ncon;
//			else
//				++ndis;
//		}
//	}
//
//
//	// concordant
//
//	if (ncon > 1)
//	{
//		const size_t   max = std::floor(static_cast<decimal_t>(ncon) * exclude_prop);
//		size_t         num = 0;
//		SortNode::List exc;
//
//		for (it = this->list.begin(); it != ti; ++it)
//		{
//			if (!(*it)->ccf[clock].good)
//				continue;
//
//			if ((*it)->sharing)
//			{
//				exc.push_back( SortNode((*it)->ccf[clock].q50, it) ); // median of CCF
//			}
//		}
//
//		exc.sort();
//		exc.reverse(); // sort highest first
//
//		for (SortNode::List::const_reference x : exc)
//		{
//			(*x.pair)->ccf[clock].good = false;
//
//			if (++num >= max)
//				break;
//		}
//	}
//
//
//	// discordant
//
//	if (ndis > 1)
//	{
//		const size_t   max = std::floor(static_cast<decimal_t>(ndis) * exclude_prop);
//		size_t         num = 0;
//		SortNode::List exc;
//
//		for (it = this->list.begin(); it != ti; ++it)
//		{
//			if (!(*it)->ccf[clock].good)
//				continue;
//
//			if (!(*it)->sharing)
//			{
//				exc.push_back( SortNode((*it)->ccf[clock].q50, it) ); // median of CCF
//			}
//		}
//
//		exc.sort(); // sort lowest first
//
//		for (SortNode::List::const_reference x : exc)
//		{
//			(*x.pair)->ccf[clock].good = false;
//
//			if (++num >= max)
//				break;
//		}
//	}
//}


// print to file

void Site::print_header(std::ostream & stream, const Param::Data param, const bool density)
{
	stream << "MarkerID Clock Filtered";

	if (density)
	{
		decimal_vector_t const & prior = param->prior; // times prior

		decimal_vector_t::const_iterator time, time_end = prior.cend();

		for (time = prior.cbegin(); time != time_end; ++time)
		{
			stream << " t" << std::fixed << std::setprecision(8) << *time;
		}
	}
	else
	{
//		stream << " Fk N_Shared N_Others PostMean PostMode PostMedian PostCI025 PostCI975 Robust Lower Upper";
		stream << " N_Concordant N_Discordant PostMean PostMode PostMedian";
	}

	stream << std::endl;
}

void Site::print(std::ostream & stream, const size_t & ne, const bool density, const bool print_adj) const
{
	if (!this->done)
		throw std::runtime_error("Invalid site result");

	for (int c = 0; c < n_clocks; ++c)
	{
		CLE const & out = (print_adj) ? this->adj[c]: this->raw[c];

		if (!out.good)
			continue;

		stream << this->focus.value << ' '; // marker id

		// clock
		if (c == 0) stream << "M ";
		if (c == 1) stream << "R ";
		if (c == 2) stream << "J ";

		if (print_adj)
			stream << "1 "; // adjusted
		else
			stream << "0 "; // raw

		if (density)
		{
			decimal_vector_t::const_iterator seq, end = out.grid.cend();

			for (seq = out.grid.cbegin(); seq != end; ++seq)
			{
				stream << ' ' << std::fixed << std::setprecision(8) << *seq;
			}
		}
		else
		{
//			stream << this->fk << ' ';

			stream << out.n_shared << ' ';
			stream << out.n_others << ' ';

			stream << std::fixed << std::setprecision(8) << out.mean   * decimal_two * static_cast<decimal_t>(ne) << ' ';
			stream << std::fixed << std::setprecision(8) << out.mode   * decimal_two * static_cast<decimal_t>(ne) << ' ';
			stream << std::fixed << std::setprecision(8) << out.median * decimal_two * static_cast<decimal_t>(ne); // << ' ';
//			stream << std::fixed << std::setprecision(8) << out.ci_025 << ' ';
//			stream << std::fixed << std::setprecision(8) << out.ci_975 << ' ';
//			stream << std::fixed << std::setprecision(8) << out.estim << ' ';
//			stream << std::fixed << std::setprecision(8) << out.lower << ' ';
//			stream << std::fixed << std::setprecision(8) << out.upper;
		}

		stream << std::endl;
	}
}



// Result queue

// construct

Queue::Hold::Hold(const size_t & _fk, const Gen::Marker::Key & _site, const Gen::Sample::Key::Vector & _share)
: fk(_fk)
, site(_site)
, share(_share)
{}

// construct for simulated results

Queue::Hold::Hold(const Gen::Marker::Key & _site)
: fk(0)
, site(_site)
, share(0)
{}


// construct

Queue::Queue(const Share::Data share, const size_t _limit, const Grid::Data _source, const Param::Data _param, const IBD::SIM::Result::Data simres, const bool random_order)
: source(_source)
, param(_param)
, limit(_limit)
, total(0)
{
	// construct for simulated results
	if (simres)
	{
		IBD::SIM::Truth::Map::const_iterator map, map_end = simres->get().cend();

		for (map = simres->get().cbegin(); map != map_end; ++map)
		{
			this->total += map->second.size();

			// hold in queue
			this->queue.emplace_back(map->first);
		}

		return;
	}


	Share::Index::Map::const_iterator map, map_end = share->get().cend();

	for (map = share->get().cbegin(); map != map_end; ++map)
	{
		const size_t fk = map->first;

		Share::Index::Sites::const_iterator idx, idx_end = map->second.sites.cend();

		for (idx = map->second.sites.cbegin(); idx != idx_end; ++idx)
		{
			const size_t n = idx->second.size();

			if (n >= Share::minimum)
			{
				//const size_t n_others = this->param->Ng - n;
				//this->total += std::min(this->param->outgroup_size, n_others);
				//const size_t n_sharer = std::min((n * (n - 1)) / 2, this->param->limit_sharers);

				this->total += std::min((n * (n - 1)) / 2, this->param->limit_sharers);
				this->total += std::min(param->outgroup_size, (this->param->Ng - n) * n);

				// hold in queue
				this->queue.emplace_back(fk, idx->first, idx->second);
			}
		}
	}

	if (random_order)
	{
		std::shuffle(this->queue.begin(), this->queue.end(), random_generator);
	}
}


// next batch

size_t Queue::next(const IBD::SIM::Result::Data simres)
{
	size_t pair_count = 0;


	// clean current batch

	if (!this->pairs.empty())
	{
		Pair::List::iterator pair = this->pairs.begin();
		Pair::List::iterator pair_end = this->pairs.end();

		while (pair != pair_end)
		{
			if ((*pair)->done)
			{
				pair = this->pairs.erase(pair);
				pair_end = this->pairs.end();
				continue;
			}

			++pair;
			++pair_count;
		}
	}

	if (!this->sites.empty())
	{
		Site::List::iterator site = this->sites.begin();
		Site::List::iterator site_end = this->sites.end();

		while (site != site_end)
		{
			if ((*site)->done)
			{
				site = this->sites.erase(site);
				site_end = this->sites.end();
				continue;
			}

			++site;
		}
	}


	// append to batch

	while (!this->queue.empty())
	{
		Hold & q = this->queue.front();

		Site::Data site;

		if (simres)
		{
			site = std::make_shared< Site >(q.site, simres); // fetch site
		}
		else
		{
			site = std::make_shared< Site >(q.fk, q.site, q.share, this->source, this->param); // make site
		}

		if (!site->done)
		{
			// copy site pairs into batch

			Pair::List::const_iterator pair, pair_end = site->list.cend();

			for (pair = site->list.cbegin(); pair != pair_end; ++pair)
			{
				(*pair)->site = site; // activate site pointer !!!

				this->pairs.push_back(*pair);

				++pair_count;
			}

			this->sites.push_back(site); // copy site into batch
		}

		this->queue.pop_front(); // remove in queue

		if (pair_count > this->limit)
		{
			break;
		}
	}


	return pair_count;
}


// total number of pairs

size_t Queue::size() const
{
	return this->total;
}


// Exection of segment detection and age inference
//

// construct

Infer::Infer(const Param::Data para, const DetectMethod detectmethod, const decimal_t max_miss, const Pair::Data target_pair, const Grid::Data grid, const HMM::Model::Data hmm_model)
: param(para)
, method(detectmethod)
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


// execute detection

void Infer::run()
{
	// lock site pointer
	
	if (this->target->site.expired())
	{
		throw std::runtime_error("Unexpected site pointer deletion");
	}
	
	Site::Data site_ptr = this->target->site.lock();
	
	
	// get data
	
	const Variant::Vector::Data a = this->source->get(this->target->pair.first.individual);
	const Variant::Vector::Data b = this->source->get(this->target->pair.second.individual);
	
	this->target->done = true;
	
	
	// check genotypes
	
	if (is_genotype<G_>(a->gen(site_ptr->focus)) || is_genotype<G_>(b->gen(site_ptr->focus)))
	{
		return;
	}
	
	//	if (this->target->sharing)
	//	{
	//		if (!is_genotype<G1>(b->gen(site_ptr->focus)))
	//			return;
	//	}
	//	else
	//	{
	//		if (!is_genotype<G0>(b->gen(site_ptr->focus)))
	//			return;
	//	}
	
	
	// determine chromosomes
	
	//	if (this->param->use_mut_clock || this->method == IBD::DETECT_FGT)
	//	{
	this->target->pair.first.chromosome = this->chr_share(a->at(site_ptr->focus));
	
	this->target->pair.second.chromosome = (this->target->sharing) ?
	this->chr_share(b->at(site_ptr->focus)):
	this->chr_other(b->at(site_ptr->focus));
	
	if (this->target->pair.first.chromosome  == CHR_VOID ||
		this->target->pair.second.chromosome == CHR_VOID)
	{
		return;
	}
	
	if (this->target->pair.first.individual  == this->target->pair.second.individual &&
		this->target->pair.first.chromosome  == this->target->pair.second.chromosome)
	{
		return;
	}
	
	//	}
	//	else
	//	{
	//		this->target->pair.first.chromosome  = UNPHASED;
	//		this->target->pair.second.chromosome = UNPHASED;
	//	}
	
	
	const hap_vector_t ha = a->hap(this->target->pair.first.chromosome);
	const hap_vector_t hb = b->hap(this->target->pair.second.chromosome);
	
	if (this->target->sharing)
	{
		if (ha.at(site_ptr->focus) != hb.at(site_ptr->focus))
			throw std::string("Unexpected allele sharing error");
	}
	else
	{
		if (ha.at(site_ptr->focus) == hb.at(site_ptr->focus))
			throw std::string("Unexpected allele non-sharing error");
	}
	
	
	// calculate missing rate
	
	this->target->missing = missing_rate(a->gen(), b->gen(), a->size());
	
	
	// detect segment
	
	this->hmm(site_ptr, ha, hb);
	
	//	if (this->target->missing < this->max_missing_rate)
	//	{
	//		switch (this->method)
	//		{
	//			case IBD::DETECT_DGT: this->dgt(site_ptr, a, b); break;
	//			case IBD::DETECT_FGT: this->fgt(site_ptr, a, b); break;
	//			case IBD::DETECT_HMM: this->hmm(site_ptr, a, b); break;
	//			case IBD::DETECT_SIM: this->sim(site_ptr, a, b); break;
	//			case IBD::DETECT_VOID: return;
	//		}
	//	}
}


// determine chromosomes

ChrType Infer::chr_share(Variant && var)
{
	if (var.is_phased())
	{
		const HapType mat = haplotype_to_index(var.hap(MATERNAL));
		const HapType pat = haplotype_to_index(var.hap(PATERNAL));

		if (is_haplotype<H1>(mat) && is_haplotype<H0>(pat))
		{
			return MATERNAL;
		}

		if (is_haplotype<H1>(pat) && is_haplotype<H0>(mat))
		{
			return PATERNAL;
		}

		if (is_haplotype<H1>(pat) && is_haplotype<H1>(mat))
		{
			return (random_coin()) ? MATERNAL: PATERNAL;
		}

		return CHR_VOID;
	}

	return UNPHASED;
}

ChrType Infer::chr_other(Variant && var)
{
	if (var.is_phased())
	{
		const HapType mat = haplotype_to_index(var.hap(MATERNAL));
		const HapType pat = haplotype_to_index(var.hap(PATERNAL));

		if (is_haplotype<H_>(mat))
		{
			if (is_haplotype<H_>(pat))
			{
				return CHR_VOID;
			}
			return PATERNAL;
		}

		if (is_haplotype<H_>(pat))
		{
			if (is_haplotype<H_>(mat))
			{
				return CHR_VOID;
			}
			return MATERNAL;
		}

		if (is_haplotype<H1>(mat))
		{
			if (is_haplotype<H1>(pat))
			{
				return CHR_VOID;
			}
			return PATERNAL;
		}

		if (is_haplotype<H1>(pat))
		{
			if (is_haplotype<H1>(mat))
			{
				return CHR_VOID;
			}
			return MATERNAL;
		}

		// choose at random
		return (random_coin()) ? MATERNAL: PATERNAL;
	}

	return UNPHASED;
}


// determine segment differences

void Infer::detect_segdiff(const Site::Data site, const hap_vector_t hap0, const hap_vector_t hap1)
{
	//	const ChrType chr0 = this->target->pair.first.chromosome;
	//	const ChrType chr1 = this->target->pair.second.chromosome;
	//
	//	if ((chr0 != MATERNAL && chr0 != PATERNAL) ||
	//		(chr1 != MATERNAL && chr1 != PATERNAL))
	//	{
	//		throw std::logic_error("Genotype data is not phased");
	//	}
	//
	//	const hap_vector_t & hap0 = a->hap(chr0);
	//	const hap_vector_t & hap1 = b->hap(chr1);
	//

	// LHS

	this->target->segdiff[LHS] = 0;

	if (site->focus.value != this->param->boundary[LHS])
	{
		for (size_t i = site->focus.value - 1, j = this->target->segment[LHS].value; i > j; --i)
		{
			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
				++this->target->segdiff[LHS];
		}

		if (this->target->segment[LHS].value == this->param->boundary[LHS])
		{
			const size_t i = this->target->segment[LHS].value;

			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
				++this->target->segdiff[LHS];
		}
	}


	// RHS

	this->target->segdiff[RHS] = 0;

	if (site->focus.value != this->param->boundary[RHS])
	{
		for (size_t i = site->focus.value + 1, j = this->target->segment[RHS].value; i < j; ++i)
		{
			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
			{
				++this->target->segdiff[RHS];
			}
		}

		if (this->target->segment[RHS].value == this->param->boundary[RHS])
		{
			const size_t i = this->target->segment[RHS].value;

			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
				++this->target->segdiff[RHS];
		}
	}
}


//void Infer::approx_segdiff(const Site::Data site, const Variant::Vector::Data a, const Variant::Vector::Data b, const bool is_shared)
//{
//	if (is_shared)
//		this->approx_segdiff_concordant(site, a, b);
//	else
//		this->detect_segdiff(site, a, b);
//	    //this->approx_segdiff_discordant(site, a, b); // not applied anymore
//}

//void Infer::approx_segdiff_concordant(const Site::Data site, const Variant::Vector::Data a, const Variant::Vector::Data b)
void Infer::approx_segdiff(const Site::Data site, const hap_vector_t hap0, const hap_vector_t hap1)
{
	const Marker::Vector & mrk = this->source->marker();

	if (mrk.at(site->focus.value).hap_count[H1] != site->fk)
		throw std::runtime_error("Inconsistent focal allele count");

	const size_t fk = site->fk;

//	const ChrType chr0 = this->target->pair.first.chromosome;
//	const ChrType chr1 = this->target->pair.second.chromosome;
//
//	if ((chr0 != MATERNAL && chr0 != PATERNAL) ||
//		(chr1 != MATERNAL && chr1 != PATERNAL))
//	{
//		throw std::logic_error("Genotype data is not phased");
//	}
//
//	const hap_vector_t & hap0 = a->hap(chr0);
//	const hap_vector_t & hap1 = b->hap(chr1);


	// LHS

	this->target->segdiff[LHS] = 0;

	if (site->focus.value != this->param->boundary[LHS])
	{
		for (size_t i = site->focus.value - 1, j = this->target->segment[LHS].value; i > j; --i)
		{
			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
			{
				if (mrk.at(i).hap_count[H1] <= fk)
					++this->target->segdiff[LHS];
			}
		}

		if (this->target->segment[LHS].value == this->param->boundary[LHS])
		{
			const size_t i = this->target->segment[LHS].value;

			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
			{
				if (mrk.at(i).hap_count[H1] <= fk)
					++this->target->segdiff[LHS];
			}
		}
	}


	// RHS

	this->target->segdiff[RHS] = 0;

	if (site->focus.value != this->param->boundary[RHS])
	{
		for (size_t i = site->focus.value + 1, j = this->target->segment[RHS].value; i < j; ++i)
		{
			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
			{
				if (mrk.at(i).hap_count[H1] <= fk)
					++this->target->segdiff[RHS];
			}
		}

		if (this->target->segment[RHS].value == this->param->boundary[RHS])
		{
			const size_t i = this->target->segment[RHS].value;

			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
			{
				if (mrk.at(i).hap_count[H1] <= fk)
					++this->target->segdiff[RHS];
			}
		}
	}
}

//void Infer::approx_segdiff_discordant(const Site::Data site, const Variant::Vector::Data a, const Variant::Vector::Data b)
//{
//	const ChrType chr0 = this->target->pair.first.chromosome;
//	const ChrType chr1 = this->target->pair.second.chromosome;
//
//	if ((chr0 != MATERNAL && chr0 != PATERNAL) ||
//		(chr1 != MATERNAL && chr1 != PATERNAL))
//	{
//		throw std::logic_error("Genotype data is not phased");
//	}
//
//	const hap_vector_t & hap0 = a->hap(chr0);
//	const hap_vector_t & hap1 = b->hap(chr1);
//
//
//	// LHS
//
//	this->target->segdiff[LHS] = 0;
//
//	if (site->focus.value != this->param->boundary[LHS])
//	{
//		for (size_t i = site->focus.value - 1, j = this->target->segment[LHS].value; i > j; --i)
//		{
//			if (is_haplotype<H1>(hap0[i]) && is_haplotype<H0>(hap1[i]))
//				++this->target->segdiff[LHS];
//		}
//
//		if (this->target->segment[LHS].value == this->param->boundary[LHS])
//		{
//			const size_t i = this->target->segment[LHS].value;
//
//			if (is_haplotype<H1>(hap0[i]) && is_haplotype<H0>(hap1[i]))
//				++this->target->segdiff[LHS];
//		}
//	}
//
//
//	// RHS
//
//	this->target->segdiff[RHS] = 0;
//
//	if (site->focus.value != this->param->boundary[RHS])
//	{
//		for (size_t i = site->focus.value + 1, j = this->target->segment[RHS].value; i < j; ++i)
//		{
//			if (hap_break(hap0[i], hap1[i]) == HAP_BREAK)
//			{
//				if (is_haplotype<H1>(hap0[i]) && is_haplotype<H0>(hap1[i]))
//					++this->target->segdiff[RHS];
//			}
//		}
//
//		if (this->target->segment[RHS].value == this->param->boundary[RHS])
//		{
//			const size_t i = this->target->segment[RHS].value;
//
//			if (is_haplotype<H1>(hap0[i]) && is_haplotype<H0>(hap1[i]))
//				++this->target->segdiff[RHS];
//		}
//	}
//}


// run methods

//void Infer::fgt(const Site::Data site, const Variant::Vector::Data a, const Variant::Vector::Data b)
//{
//	// ibd detection
//
//	site->frequency(this->source); // calculate sharer subset allele frequency
//
//	if (this->target->sharing)
//	{
//		FGT::Algorithm algorithm(a, b);
//
//		this->target->segment = algorithm.detect(site->focus);
//	}
//	else
//	{
//		FGT::Algorithm algorithm(a, this->target->pair.first.chromosome,
//								 b, this->target->pair.second.chromosome,
//								 site->freq);
//
//		this->target->segment = algorithm.detect(site->focus);
//	}
//
//	// age estimation
//
//	Density age(this->param, this->target->sharing, site->focus, this->target->segment);
//
//	//if (this->param->use_mut_clock)
//	//{
//	if (this->param->use_tree_consistency)
//		this->approx_segdiff(site, a, b, this->target->sharing);
//	else
//		this->detect_segdiff(site, a, b);
//
//	age.differences(this->target->segdiff);
//	//}
//
//	if (this->param->run_mut_clock) this->target->ccf[MUT_CLOCK] = age.estimate(MUT_CLOCK);
//	if (this->param->run_rec_clock) this->target->ccf[REC_CLOCK] = age.estimate(REC_CLOCK);
//	if (this->param->run_cmb_clock) this->target->ccf[CMB_CLOCK] = age.estimate(CMB_CLOCK);
//}
//
//void Infer::dgt(const Site::Data site, const Variant::Vector::Data a, const Variant::Vector::Data b)
//{
//	// ibd detection
//
//	DGT::Algorithm algorithm(a, b);
//
//	this->target->segment = algorithm.detect(site->focus);
//
//
//	// age estimation
//
//	Density age(this->param, this->target->sharing, site->focus, this->target->segment);
//
//	//if (this->param->use_mut_clock)
//	//{
//	if (this->param->use_tree_consistency)
//		this->approx_segdiff(site, a, b, this->target->sharing);
//	else
//		this->detect_segdiff(site, a, b);
//
//	age.differences(this->target->segdiff);
//	//}
//
//	if (this->param->run_mut_clock) this->target->ccf[MUT_CLOCK] = age.estimate(MUT_CLOCK);
//	if (this->param->run_rec_clock) this->target->ccf[REC_CLOCK] = age.estimate(REC_CLOCK);
//	if (this->param->run_cmb_clock) this->target->ccf[CMB_CLOCK] = age.estimate(CMB_CLOCK);
//}

void Infer::hmm(const Site::Data site, const hap_vector_t a, const hap_vector_t b)
{
	// ibd detection

	HMM::Algorithm algorithm(a, b, this->model, !this->target->sharing);

	//	if (this->model->do_iterative)
	//		this->target->segment = algorithm.detect_iterative(site->fk, site->focus);
	//	else
	this->target->segment = algorithm.detect(site->fk, site->focus);


	// age estimation

	Density age(this->param, this->target->sharing, site->focus, this->target->segment);

	
	if (this->param->use_tree_consistency && this->target->sharing)
		this->approx_segdiff(site, a, b);
	else
		this->detect_segdiff(site, a, b);

	
	age.differences(this->target->segdiff);
	

	if (this->param->use_post_prob)
	{
		age.probability(algorithm.posterior(LHS, HMM::NON_STATE, true), algorithm.posterior(RHS, HMM::NON_STATE, true));
	}

	if (this->param->run_mut_clock) this->target->ccf[MUT_CLOCK] = age.estimate(MUT_CLOCK);
	if (this->param->run_rec_clock) this->target->ccf[REC_CLOCK] = age.estimate(REC_CLOCK);
	if (this->param->run_cmb_clock) this->target->ccf[CMB_CLOCK] = age.estimate(CMB_CLOCK);
}

//void Infer::sim(const Site::Data site, const Variant::Vector::Data a, const Variant::Vector::Data b)
//{
//	// age estimation
//
//	Density age(this->param, this->target->sharing, site->focus, this->target->segment);
//
//	//if (this->param->use_mut_clock)
//	//{
//	if (this->param->use_tree_consistency)
//		this->approx_segdiff(site, a, b, this->target->sharing);
//	else
//		this->detect_segdiff(site, a, b);
//
//	age.differences(this->target->segdiff);
//	//}
//
//	if (this->param->run_mut_clock) this->target->ccf[MUT_CLOCK] = age.estimate(MUT_CLOCK);
//	if (this->param->run_rec_clock) this->target->ccf[REC_CLOCK] = age.estimate(REC_CLOCK);
//	if (this->param->run_cmb_clock) this->target->ccf[CMB_CLOCK] = age.estimate(CMB_CLOCK);
//}

