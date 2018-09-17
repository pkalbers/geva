//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "GenShare.hpp"


using namespace Gen;


// constant char sequence as identifier for binary file navigation
const char Share::checkpoint[4] = { 0x25, 0x50, 0x4b, 0x41 };


// construct

Share::Share()
: region(0, 0)
, max_sites(0)
, max_pairs(0)
, good(false)
{}


// create new index for target genotype count

bool Share::detect(const target_t & target_fk, const Grid::Data grid)
{
	if (this->good)
	{
		throw std::runtime_error("Sharing index already generated");
	}
	
	
	// limit region
	const bool in_region = (this->region.first != 0 || this->region.second != 0);
	
	
	target_t::const_iterator fk, fk_end = target_fk.cend();
	
	for (fk = target_fk.cbegin(); fk != fk_end; ++fk)
	{
		const size_t target = *fk;
		
		if (target < minimum)
			continue;
		
		
		Index index;
		
		index.fk = target;
		
		
		// detect markers
		
		size_t count = 0;
		
		Marker::Iterator marker, marker_end = grid->marker().cend();
		
		for (marker = grid->marker().cbegin(); marker != marker_end; ++marker)
		{
			if (in_region)
			{
				const size_t pos = marker->position;
				
				if (pos < this->region.first)
					continue;
				
				if (pos > this->region.second)
					break;
			}
			
			if (marker->hap_count[H1] == target && marker->gen_count[G1] == target) // !!! condition: allele and genotype count
			{
				index.sites[ marker->index ] = Sample::Key::Vector();
				index.sites[ marker->index ].reserve(target);
				
				++count;
			}
		}
		
		
		if (count > 0)
		{
			// limit sites
			if (this->max_sites > 0)
			{
				this->subset_sites_index(index, this->max_sites);
			}
			
			// move index into table
			const std::pair< Index::Map::iterator, bool > out = this->table.emplace(target, std::move(index));
			
			if (!out.second)
				throw std::runtime_error("Unable to store index in memory");
		}
	}
	
	
	this->detect_share(grid);
	
	this->create_pairs_table();
	
	
	// limit pairs
	if (this->max_pairs > 0)
	{
		this->subset_pairs_table();
		this->remake_sites_table();
	}
	
	
	this->clean_table();
	
	this->good = true;
	
	return (!this->table.empty());
}


// create index for input list of target sites

size_t Share::select(const std::set<size_t> & position, const Grid::Data grid)
{
	if (this->good)
	{
		throw std::runtime_error("Sharing index already generated");
	}
	
	
	size_t count = 0;
	
	Marker::Iterator marker = grid->marker().cbegin(), marker_end = grid->marker().cend();
	
	
	// limit region
	const bool in_region = (this->region.first != 0 || this->region.second != 0);
	
	
	// walkabout focal positions
	
	std::set<size_t>::const_iterator pos, end = position.cend();
	
	for (pos = position.cbegin(); pos != end; ++pos)
	{
		const size_t focus = *pos;
		
		if (in_region)
		{
			if (focus < this->region.first)
				continue;
			
			if (focus > this->region.second)
				break;
		}
		
		// walkabout markers
		for (; marker != marker_end; ++marker)
		{
			if (marker->position == focus)
			{
				const size_t target = marker->hap_count[H1]; // allele count at focal site
				
				if (target >= minimum)
				{
					// find/make corresponding index
					Index & index = this->table[ target ];
					
					index.fk = target;
					
					// insert focal site
					index.sites[ marker->index ] = Sample::Key::Vector();
					index.sites[ marker->index ].reserve( target );
					
					++count;
				}
				
				break;
			}
		}
	}
	
//	printf("\n %lu \n", count);
	
	if (count > 0)
	{
		count = 0;
		
		
		// limit sites
		if (this->max_sites > 0)
		{
			this->subset_sites_table();
		}
		
		
		this->detect_share(grid);
		
		
		// walkabout share targets
		
		Index::Map::iterator index, index_end = this->table.end();
		
		
		for (index = this->table.begin(); index != index_end; ++index)
		{
			count += index->second.sites.size();
		}
//		printf("\n %lu \n", count);
		count = 0;
		
		
		for (index = this->table.begin(); index != index_end; ++index)
		{
			create_pairs_index(index->second);
			
			// limit pairs
			if (this->max_pairs > 0)
			{
				subset_pairs_index(index->second, this->max_pairs);
				remake_sites_index(index->second);
			}
			
			count += index->second.sites.size();
		}
	}
	
//	printf("\n %lu \n", count);
	
	this->clean_table();
	
	this->good = true;
	
	return count;
}


// access table

Share::Index::Map const & Share::get() const
{
	if (!this->good)
	{
		throw std::runtime_error("Sharing index not generated");
	}
	
	return this->table;
}


// print size of index to stream

void Share::print(const size_t & N, std::ostream & stream) const
{
	const double n = static_cast<double>(N);
	const int    p = static_cast<int>(std::ceil(std::log10(n / 100.0)));
	
	size_t total_sites = 0;
	size_t total_pairs = 0;
	
	
	// header
	stream << std::right << std::setw(10) << "# shared" << " ";
	stream << std::right << std::setw(10) << "% freq."  << " ";
	stream << std::right << std::setw(10) << "# sites"  << " ";
	stream << std::right << std::setw(10) << "# pairs"  << std::endl;
	
	
	Index::Map::const_iterator index, index_end = this->table.cend();
	
	for (index = this->table.cbegin(); index != index_end; ++index)
	{
		const size_t n_sites = index->second.sites.size();
		const size_t n_pairs = index->second.pairs.size();
		
		stream << std::fixed;
		stream << std::right << std::setw(10) << index->first << " ";
		stream << std::right << std::setw(10) << std::setprecision(std::min(8, std::max(p, 1))) << ((index->first / n) * 100.0) << " ";
		stream << std::right << std::setw(10) << n_sites << " ";
		stream << std::right << std::setw(10) << n_pairs << std::endl;
		
		total_sites += n_sites;
		total_pairs += n_pairs;
	}
	
	
	// footer
	stream << std::string(10, ' ') << " " << std::string(10, ' ') << " ";
	stream << std::string(10, '-') << " " << std::string(10, '-') << std::endl;
	stream << " " << std::left << std::setw(20) << "Total:" << " ";
	stream << std::right << std::setw(10) << total_sites << " ";
	stream << std::right << std::setw(10) << total_pairs << std::endl;
}


// detect share across samples

void Share::detect_share(const Grid::Data grid)
{
	Sample::Iterator sample, sample_end = grid->sample().cend();
	
	// walkabout samples
	for (sample = grid->sample().cbegin(); sample != sample_end; ++sample)
	{
		Grid::Vector data = grid->read(sample->index);
		
		
		// walkabout indices
		
		Index::Map::iterator index, index_end = this->table.end();
		
		for (index = this->table.begin(); index != index_end; ++index)
		{
			const size_t target = index->first;
			
			
			// walkabout sites
			
			Index::Sites::iterator site, site_end = index->second.sites.end();
			
			for (site = index->second.sites.begin(); site != site_end; ++site)
			{
				//if (is_genotype<G1>(data[ site->first.value ])) // !!! condition: heterozygous target site
				if (is_genotype<G1>(data[ site->first.value ]))
				{
					if (site->second.size() == target)
					{
						throw std::logic_error("Unexpected number of individuals");
					}
					
					site->second.push_back(sample->index);
				}
				else if (is_genotype<G2>(data[ site->first.value ]))
				{
					if (site->second.size() == target)
					{
						throw std::logic_error("Unexpected number of individuals");
					}
					
					site->second.push_back(sample->index);
					site->second.push_back(sample->index);
				}
			}
		}
	}
}


// make sharer pairs

void Share::create_pairs_index(Index & index)
{
	index.pairs.clear();
	
	Index::Sites::iterator site = index.sites.begin(), site_end = index.sites.end();
	
	while (site != site_end)
	{
		const size_t n = site->second.size();
		
		if (n < minimum)
		{
			site = index.sites.erase(site);
			continue;
		}
		
		// make each pair
		for (size_t i = 0; i < n - 1; ++i)
		{
			for (size_t j = i + 1; j < n; ++j)
			{
				Sample::Key::Pair pair(site->second[i], site->second[j]);
				
				index.pairs[ pair ].insert(site->first);
			}
		}
		
		++site;
	}
}

void Share::create_pairs_table()
{
	Index::Map::iterator index, index_end = this->table.end();
	
	for (index = this->table.begin(); index != index_end; ++index)
	{
		create_pairs_index(index->second);
	}
}


// re-generate sites from pairs

void Share::remake_sites_index(Index & index)
{
	index.sites.clear();
	
	Index::Pairs::const_iterator pair, pair_end = index.pairs.cend();
	
	for (pair = index.pairs.cbegin(); pair != pair_end; ++pair)
	{
		const Sample::Key a = pair->first.first;
		const Sample::Key b = pair->first.second;
		
		Marker::Key::Set::const_iterator site, site_end = pair->second.cend();
		
		for (site = pair->second.cbegin(); site != site_end; ++site)
		{
			const Marker::Key key = *site;
			
			if (index.sites.count(key) == 0)
			{
				index.sites[ key ] = Sample::Key::Vector();
				index.sites[ key ].reserve(index.fk);
			}
			
			index.sites[ key ].push_back(a);
			index.sites[ key ].push_back(b);
		}
	}
	
	
	// sort sharers per site
	
	Index::Sites::iterator site, site_end = index.sites.end();
	
	for (site = index.sites.begin(); site != site_end; ++site)
	{
		std::sort(site->second.begin(), site->second.end());
	}
}

void Share::remake_sites_table()
{
	Index::Map::iterator index, index_end = this->table.end();
	
	for (index = this->table.begin(); index != index_end; ++index)
	{
		remake_sites_index(index->second);
	}
}


// random subset of sites or pairs

void Share::subset_sites_index(Index & index, const size_t & max)
{
	const size_t len = index.sites.size();
	
	if (len <= max)
	{
		return;
	}
	
	
	std::vector< Index::Sites::iterator > its;
	its.reserve(len);
	
	Index::Sites::iterator it, ti = index.sites.end();
	
	for (it = index.sites.begin(); it != ti; ++it)
	{
		its.push_back(it);
	}
	
	std::shuffle(its.begin(), its.end(), random_generator);
	
	
	Index::Sites sub_sites;
	
	for (size_t i = 0; i < max; ++i)
	{
		sub_sites[ its[i]->first ].swap(its[i]->second);
	}
	
	index.sites.swap(sub_sites);
}

void Share::subset_sites_table()
{
	Index::Map::iterator index, index_end = this->table.end();
	
	for (index = this->table.begin(); index != index_end; ++index)
	{
		subset_sites_index(index->second, this->max_sites);
	}
}


// random subset of pairs

void Share::subset_pairs_index(Index & index, const size_t & max)
{
	const size_t len = index.pairs.size();
	
	if (len <= max)
	{
		return;
	}
	
	std::vector< Index::Pairs::iterator > its;
	its.reserve(len);
	
	Index::Pairs::iterator it, ti = index.pairs.end();
	
	for (it = index.pairs.begin(); it != ti; ++it)
	{
		its.push_back(it);
	}
	
	std::shuffle(its.begin(), its.end(), random_generator);
	
	
	Index::Pairs sub_pairs;
	
	for (size_t i = 0; i < max; ++i)
	{
		sub_pairs[ its[i]->first ].swap(its[i]->second);
	}
	
	index.pairs.swap(sub_pairs);
}

void Share::subset_pairs_table()
{
	Index::Map::iterator index, index_end = this->table.end();
	
	for (index = this->table.begin(); index != index_end; ++index)
	{
		this->subset_pairs_index(index->second, this->max_pairs);
	}
}


// remove empty indices in table

void Share::clean_table()
{
	Index::Map::iterator index = this->table.begin(), index_end = this->table.end();
	
	while (index != index_end)
	{
		if (index->second.sites.size() == 0 || index->second.pairs.size() == 0)
		{
			index = this->table.erase(index);
			continue;
		}
		
		++index;
	}
}



//
// load and save functions
//

// load index map from binary file

bool Share::load(const std::string & filename, const Grid::Data grid)
{
	if (this->good)
	{
		throw std::runtime_error("Sharing index already generated");
	}
	
	
	// limit region
	const bool in_region = (this->region.first != 0 || this->region.second != 0);
	
	
	Binary bin(filename, Binary::mode::READ);
	
	const size_t n_index = this->load_header(bin, grid);
	
	for (size_t i = 0; i < n_index; ++i)
	{
		bin.match<char>(checkpoint, 4);
		
		Index index;
		
		index.fk = bin.read<Bin4, size_t>(); // fk
		
		const size_t n_sites = bin.read<Bin4, size_t>(); // number of sites
		const size_t n_pairs = bin.read<Bin4, size_t>(); // number of pairs
		
		
		// sites
		
		bin.match<char>(checkpoint, 4);
		
		for (size_t k = 0; k < n_sites; ++k)
		{
			index.sites.insert(this->load_index_site(bin));
		}
		
		
		// pairs
		
		bin.match<char>(checkpoint, 4);
		
		for (size_t k = 0; k < n_pairs; ++k)
		{
			index.pairs.insert(this->load_index_pair(bin));
		}
		
		
		bin.match<char>(checkpoint, 4);
		
		
		if (in_region)
		{
			size_t count = 0;
			
			Index::Sites::iterator site = index.sites.begin(), site_end = index.sites.end();
			
			while (site != site_end)
			{
				const size_t pos = grid->marker(site->first).position;
				
				if (pos < this->region.first ||
					pos > this->region.second)
				{
					site = index.sites.erase(site);
					++count;
					continue;
				}
				
				++site;
			}
			
			if (count > 0)
			{
				create_pairs_index(index);
			}
		}
		
		
		// limit sites
		if (this->max_sites > 0)
		{
			subset_sites_index(index, this->max_sites);
			create_pairs_index(index);
		}
		
		
		// limit pairs
		if (this->max_pairs > 0)
		{
			subset_pairs_index(index, this->max_pairs);
			remake_sites_index(index);
		}
		
		
		// move index into table
		if (index.sites.size() > 0)
		{
			const std::pair< Index::Map::iterator, bool > out = this->table.emplace(index.fk, std::move(index));
			
			if (!out.second)
				throw std::runtime_error("Unable to store index in memory");
		}
	}
	
	bin.match<char>(checkpoint, 4);
	
	
	this->clean_table();
	
	this->good = true;
	
	return (!this->table.empty());
}


size_t Share::load_header(Binary & bin, const Grid::Data grid)
{
	bin.match<char>(checkpoint, 4);
	
	const size_t sample_size = bin.read<Bin4, size_t>(); // number of samples
	const size_t marker_size = bin.read<Bin4, size_t>(); // number of markers
	
	if (sample_size != grid->sample_size())
		throw std::invalid_argument("Sample size does not match with data");
	
	if (marker_size != grid->marker_size())
		throw std::invalid_argument("Number of variants does not match with data");
	
	return bin.read<Bin4, size_t>(); // number of fk indices
}

Share::Index::Sites::value_type Share::load_index_site(Binary & bin) const
{
	bin.match<char>(checkpoint, 4);
	
	const size_t site = bin.read<Bin4, size_t>(); // focal site
	
	const size_t n_sharers = bin.read<Bin4, size_t>(); // number of sharers
	
	Sample::Key::Vector sharers(n_sharers);
	
	for (size_t i = 0; i < n_sharers; ++i)
	{
		sharers[i] = bin.read<Bin4, size_t>(); // sharer
	}
	
	return Share::Index::Sites::value_type(site, std::move(sharers));
}

Share::Index::Pairs::value_type Share::load_index_pair(Binary & bin) const
{
	bin.match<char>(checkpoint, 4);
	
	Sample::Key::Pair pair;
	
	pair.first  = bin.read<Bin4, size_t>(); // first sharer
	pair.second = bin.read<Bin4, size_t>(); // second sharer
	
	const size_t n_sites = bin.read<Bin4, size_t>(); // number of focal sites
	
	Marker::Key::Set sites;
	
	for (size_t i = 0; i < n_sites; ++i)
	{
		sites.emplace(bin.read<Bin4, size_t>()); // focal site
	}
	
	return Index::Pairs::value_type(std::move(pair), std::move(sites));
}


// save index map to binary file

void Share::save(const std::string & filename, const Grid::Data grid) const
{
	if (!this->good)
	{
		throw std::runtime_error("Sharing index was not generated");
	}
	
	Binary bin(filename, Binary::mode::WRITE);
	
	this->save_header(bin, grid);
	
	Index::Map::const_iterator map, end = this->table.cend();
	
	for (map = this->table.cbegin(); map != end; ++map)
	{
		bin.write<char>(checkpoint, 4);
		bin.write<Bin4>(map->first); // fk
		bin.write<Bin4>(map->second.sites.size()); // number of sites
		bin.write<Bin4>(map->second.pairs.size()); // number of pairs
		
		
		// sites
		
		bin.write<char>(checkpoint, 4);
		
		Index::Sites::const_iterator site, site_end = map->second.sites.cend();
		
		for (site = map->second.sites.cbegin(); site != site_end; ++site)
		{
			this->save_index_site(bin, site);
		}
		
		
		// pairs
		
		bin.write<char>(checkpoint, 4);
		
		Index::Pairs::const_iterator pair, pair_end = map->second.pairs.cend();
		
		for (pair = map->second.pairs.cbegin(); pair != pair_end; ++pair)
		{
			this->save_index_pair(bin, pair);
		}
		
		
		bin.write<char>(checkpoint, 4);
	}
	
	bin.write<char>(checkpoint, 4);
}


void Share::save_header(Binary & bin, const Grid::Data grid) const
{
	bin.write<char>(checkpoint, 4);
	bin.write<Bin4>(grid->sample_size()); // number of samples
	bin.write<Bin4>(grid->marker_size()); // number of markers
	bin.write<Bin4>(this->table.size()); // number of fk indices
}

void Share::save_index_site(Binary & bin, const Index::Sites::const_iterator & site) const
{
	bin.write<char>(checkpoint, 4);
	bin.write<Bin4>(site->first.value); // focal site
	bin.write<Bin4>(site->second.size()); // number of sharers
	
	Sample::Key::Vector::const_iterator it, ti = site->second.cend();
	
	for (it = site->second.cbegin(); it != ti; ++ it)
	{
		bin.write<Bin4>(it->value); // sharer
	}
}

void Share::save_index_pair(Binary & bin, const Index::Pairs::const_iterator & pair) const
{
	bin.write<char>(checkpoint, 4);
	bin.write<Bin4>(pair->first.first.value); // first sharer
	bin.write<Bin4>(pair->first.second.value); // second sharer
	bin.write<Bin4>(pair->second.size()); // number of focal sites
	
	Marker::Key::Set::const_iterator it, ti = pair->second.cend();
	
	for (it = pair->second.cbegin(); it != ti; ++ it)
	{
		bin.write<Bin4>(it->value); // focal site
	}
}

