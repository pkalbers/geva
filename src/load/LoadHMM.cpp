//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "LoadHMM.hpp"


using namespace IBD;
using namespace HMM;
using namespace Gen;


// construct
LoadHMM::LoadHMM(const size_t & effective_size, const size_t & sample_size)
: Ne(effective_size)
, Nh(sample_size)
, inits_good(false)
, emiss_good(false)
, dists_good(false)
, done(false)
{}


// read initial state probabilties from file and generate matrix

void LoadHMM::make_initial(const Grid::Data grid, const std::string & filename, std::ofstream * file_ptr)
{
	static constexpr decimal_t neg_one = decimal_one - 1e-03;
	static constexpr decimal_t pos_one = decimal_one + 1e-03;
	
	static const size_t header_fields = 5;
	static const std::string header[header_fields] = {
		"Frequency",
		"CON_NON",
		"CON_IBD",
		"DIS_NON",
		"DIS_IBD"};
	
	
	if (this->inits_good)
		throw std::runtime_error("Initial state probabilities already generated");
	
	
	Reader stream(filename);
	
	// read header
	if (stream.next())
	{
		Reader::Current line = stream.line();
		
		while (line.split())
		{
			Reader::Current field = line.field();
			
			if (field.number >= header_fields)
			{
				throw std::invalid_argument("Invalid file format [HMM, initial probabilities]\n" + filename +
											"\n Header contains more than " + std::to_string(header_fields) + " fields");
			}
			
			if (field.convert<std::string>() != header[field.number])
			{
				throw std::invalid_argument("Invalid file format [HMM, initial probabilities]\n" + filename +
											"\n Header requires '" + header[field.number] + "' in field " + std::to_string(field.number + 1));
			}
		}
	}
	
	
	// get sample size of current dataset
	const size_t N = grid->sample().size() * 2; // number of sampled haplotypes
	
	if (N != this->Nh)
	{
		throw std::runtime_error("Unexpected sample size detected");
	}
	
	
	std::set<std::string> unique;
	
	// read content
	while (stream.next())
	{
		Reader::Current line = stream.line();
		
		decimal_t frq = decimal_nil;
		
		decimal_t con_non = decimal_nil;
		decimal_t con_ibd = decimal_nil;
		decimal_t con_sum = decimal_nil;
		
		decimal_t dis_non = decimal_nil;
		decimal_t dis_ibd = decimal_nil;
		decimal_t dis_sum = decimal_nil;
		
		int i = 0;
		
		while (line.split())
		{
			Reader::Current field = line.field();
			
			switch (field.number)
			{
				case 0: // Frequency
				{
					frq = field.convert<decimal_t>();
					
					// check duplicate
					const size_t size = unique.size();
					unique.insert(field.convert<std::string>());
					if (size == unique.size())
					{
						throw std::invalid_argument("Duplicate frequency, line " + std::to_string(line.number + 1));
					}
					
					break;
				}
				case 1: // CON_NON
				{
					con_non = decimal_err + field.convert<decimal_t>();
					con_sum += con_non;
					break;
				}
				case 2: // CON_IBD
				{
					con_ibd = decimal_err + field.convert<decimal_t>();
					con_sum += con_ibd;
					break;
				}
				case 3: // DIS_NON
				{
					dis_non = decimal_err + field.convert<decimal_t>();
					dis_sum += dis_non;
					break;
				}
				case 4: // DIS_IBD
				{
					dis_ibd = decimal_err + field.convert<decimal_t>();
					dis_sum += dis_ibd;
					break;
				}
				
				default:
				{
					throw std::invalid_argument("More fields that expected, line " + std::to_string(line.number + 1));
					break;
				}
			}
			
			++i;
		}
		
		if (i != header_fields)
		{
			throw std::invalid_argument("Fewer fields that expected, line " + std::to_string(line.number + 1));
		}
		
		// check sum of probabilities
		if (con_sum < neg_one || con_sum > pos_one)
		{
			throw std::invalid_argument("Concordant probabilities do not sum to 1, line " + std::to_string(line.number + 1));
		}
		if (dis_sum < neg_one || dis_sum > pos_one)
		{
			throw std::invalid_argument("Discordant probabilities do not sum to 1, line " + std::to_string(line.number + 1));
		}
		
		
		const size_t num = lround(frq * static_cast<decimal_t>(this->Nh));
		
		// store input, ensure correct proportions
		Model::inits_type point_con, point_dis;
		
		point_con[ NON_STATE ] = con_non / con_sum;
		point_con[ IBD_STATE ] = con_ibd / con_sum;
		
		point_dis[ NON_STATE ] = dis_non / dis_sum;
		point_dis[ IBD_STATE ] = dis_ibd / dis_sum;
		
		this->input_inits_con[ num ] = point_con;
		this->input_inits_dis[ num ] = point_dis;
	}
	
	
	// insert extremes
	
	Model::inits_type extreme0, extreme1;
	
	extreme0[ NON_STATE ] = 1;
	extreme0[ IBD_STATE ] = 0;
	
	extreme1[ NON_STATE ] = 0;
	extreme1[ IBD_STATE ] = 1;
	
	this->input_inits_con[0] = extreme0;
	this->input_inits_con[N] = extreme1;
	
	extreme0[ NON_STATE ] = 1;
	extreme0[ IBD_STATE ] = 0;
	
	extreme1[ NON_STATE ] = 1;
	extreme1[ IBD_STATE ] = 0;
	
	this->input_inits_dis[0] = extreme0;
	this->input_inits_dis[N] = extreme1;
	
	
	// approximate for each genotype count
	
	// concordant
	for (size_t k = 1; k < N; ++k)
	{
		inits_input::iterator i1 = this->input_inits_con.lower_bound(k);
		
		if (i1->first == k)
		{
			continue;
		}
		
		inits_input::iterator i0 = std::prev(i1);
		
		decimal_t ibd = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ IBD_STATE ], i1->second[ IBD_STATE ], decimal_nil, decimal_one);
		decimal_t non = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ NON_STATE ], i1->second[ NON_STATE ], decimal_nil, decimal_one);
		
		decimal_t sum = ibd + non;
		
		Model::inits_type point;
		
		point[IBD_STATE] = ibd / sum;
		point[NON_STATE] = non / sum;
		
		this->input_inits_con[ k ] = point;
	}
	
	// discordant
	for (size_t k = 1; k < N; ++k)
	{
		inits_input::iterator i1 = this->input_inits_dis.lower_bound(k);
		
		if (i1->first == k)
		{
			continue;
		}
		
		inits_input::iterator i0 = std::prev(i1);
		
		decimal_t ibd = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ IBD_STATE ], i1->second[ IBD_STATE ], decimal_nil, decimal_one);
		decimal_t non = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ NON_STATE ], i1->second[ NON_STATE ], decimal_nil, decimal_one);
		
		decimal_t sum = ibd + non;
		
		Model::inits_type point;
		
		point[IBD_STATE] = ibd / sum;
		point[NON_STATE] = non / sum;
		
		this->input_inits_dis[ k ] = point;
	}
	
	
	// fill in for each site
	
	this->inits_con.resize(grid->marker().size());
	this->inits_dis.resize(grid->marker().size());
	
	Marker::Vector::const_iterator M, M_end = grid->marker().cend();
	
	for (M = grid->marker().cbegin(); M != M_end; ++M)
	{
		this->inits_con[ M->index.value ] = this->input_inits_con[ M->hap_count[H1] ];
		this->inits_dis[ M->index.value ] = this->input_inits_dis[ M->hap_count[H1] ];
	}
	
	// print
	if (file_ptr != nullptr)
	{
		this->print_initial(*file_ptr);
	}
	
	// finish
	this->inits_good = true;
}


// expected initial probabilities

void LoadHMM::make_initial(const Gen::Grid::Data grid)
{
	// get sample size of current dataset
	const size_t N = grid->sample().size() * 2; // number of sampled haplotypes
	
	for (size_t k = 1; k < N; ++k)
	{
		Model::inits_type point;
		
		point[IBD_STATE] = 1;
		point[NON_STATE] = 0;
		
		this->input_inits_con[ k ] = point;
		this->input_inits_dis[ k ] = point;
	}
	
	// fill in for each site
	
	this->inits_con.resize(grid->marker().size());
	this->inits_dis.resize(grid->marker().size());
	
	Marker::Vector::const_iterator M, M_end = grid->marker().cend();
	
	for (M = grid->marker().cbegin(); M != M_end; ++M)
	{
		this->inits_con[ M->index.value ] = this->input_inits_con[ M->hap_count[H1] ];
		this->inits_dis[ M->index.value ] = this->input_inits_dis[ M->hap_count[H1] ];
	}
	
	// finish
	this->inits_good = true;
}


// read emission probabilties from file and generate matrix

void LoadHMM::make_emission(const Grid::Data grid, const std::string & filename, std::ofstream * file_ptr)
{
	static const size_t header_fields = 7;
	static const std::string header[header_fields] = {
		"Frequency",
		"NON_00",
		"NON_01",
		"NON_11",
		"IBD_00",
		"IBD_01",
		"IBD_11"};
	
	
	if (this->emiss_good)
		throw std::runtime_error("Emission probabilities already generated");
	
	
	Reader stream(filename);
	
	// read header
	if (stream.next())
	{
		Reader::Current line = stream.line();
		
		while (line.split())
		{
			Reader::Current field = line.field();
			
			if (field.number >= header_fields)
			{
				throw std::invalid_argument("Invalid file format [HMM, emission probabilities]\n" + filename +
											"\n Header contains more than " + std::to_string(header_fields) + " fields");
			}
			
			if (field.convert<std::string>() != header[field.number])
			{
				throw std::invalid_argument("Invalid file format [HMM, emission probabilities]\n" + filename +
											"\n Header requires '" + header[field.number] + "' in field " + std::to_string(field.number + 1));
			}
		}
	}
	
	
	// get sample size of current dataset
	const size_t N = grid->sample().size() * 2; // number of sampled haplotypes
	
	if (N != this->Nh)
	{
		throw std::runtime_error("Unexpected sample size detected");
	}
	
	
	std::set<std::string> unique;
	
	// read content
	while (stream.next())
	{
		Reader::Current line = stream.line();
		
		decimal_t frq = 0.0;
		
		decimal_t non_00 = 0.0;
		decimal_t non_01 = 0.0;
		decimal_t non_11 = 0.0;
		
		decimal_t ibd_00 = 0.0;
		decimal_t ibd_01 = 0.0;
		decimal_t ibd_11 = 0.0;
		
		decimal_t non_sum = 0.0;
		decimal_t ibd_sum = 0.0;
		
		
		int i = 0;
		
		while (line.split())
		{
			Reader::Current field = line.field();
			
			switch (field.number)
			{
				case 0: // Frequency
				{
					frq = field.convert<decimal_t>();
					
					// check duplicate
					const size_t size = unique.size();
					unique.insert(field.convert<std::string>());
					if (size == unique.size())
					{
						throw std::invalid_argument("Duplicate frequency, line " + std::to_string(line.number + 1));
					}
					
					break;
				}
					
				case 1: { non_00 = decimal_err + field.convert<decimal_t>(); non_sum += non_00; break; }
				case 2: { non_01 = decimal_err + field.convert<decimal_t>(); non_sum += non_01; break; }
				case 3: { non_11 = decimal_err + field.convert<decimal_t>(); non_sum += non_11; break; }
					
				case 4: { ibd_00 = decimal_err + field.convert<decimal_t>(); ibd_sum += ibd_00; break; }
				case 5: { ibd_01 = decimal_err + field.convert<decimal_t>(); ibd_sum += ibd_01; break; }
				case 6: { ibd_11 = decimal_err + field.convert<decimal_t>(); ibd_sum += ibd_11; break; }
				
				default:
				{
					throw std::invalid_argument("More fields that expected, line " + std::to_string(line.number + 1));
					break;
				}
			}
			
			++i;
		}
		
		if (i != header_fields)
		{
			throw std::invalid_argument("Fewer fields that expected, line " + std::to_string(line.number + 1));
		}
		
		// check sum of probabilities
		if (non_sum < 0.999 || non_sum > 1.001)
		{
			throw std::invalid_argument("State probabilities for NON state do not sum to 1, line " + std::to_string(line.number + 1));
		}
		if (ibd_sum < 0.999 || ibd_sum > 1.001)
		{
			throw std::invalid_argument("State probabilities for IBD state do not sum to 1, line " + std::to_string(line.number + 1));
		}
		
		
		const size_t num = lround(frq * this->Nh);
		
		// store input, ensure correct proportions
		Model::emiss_type point;
		
		point[ NON_STATE ][ H00 ] = non_00 / non_sum;
		point[ NON_STATE ][ H01 ] = non_01 / non_sum;
		point[ NON_STATE ][ H11 ] = non_11 / non_sum;
		
		point[ IBD_STATE ][ H00 ] = ibd_00 / ibd_sum;
		point[ IBD_STATE ][ H01 ] = ibd_01 / ibd_sum;
		point[ IBD_STATE ][ H11 ] = ibd_11 / ibd_sum;
		
		this->input_emiss[ num ] = point;
	}
	
	
	// insert extremes
	
	Model::emiss_type extreme0, extreme1;
	
	extreme0[ NON_STATE ] = { 1, 0, 0 };
	extreme0[ IBD_STATE ] = { 1, 0, 0 };
	
	extreme1[ NON_STATE ] = { 0, 0, 1 };
	extreme1[ IBD_STATE ] = { 0, 0, 1 };
	
	this->input_emiss[0] = extreme0;
	this->input_emiss[N] = extreme1;
	
	
	// approximate for each genotype count
	
	for (size_t k = 1; k < N; ++k)
	{
		emiss_input::iterator i1 = this->input_emiss.lower_bound(k);
		
		if (i1->first == k)
		{
			continue;
		}
		
		emiss_input::iterator i0 = std::prev(i1);
		
		const decimal_t non_00 = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ NON_STATE ][ H00 ], i1->second[ NON_STATE ][ H00 ], decimal_nil, decimal_one);
		const decimal_t non_01 = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ NON_STATE ][ H01 ], i1->second[ NON_STATE ][ H01 ], decimal_nil, decimal_one);
		const decimal_t non_11 = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ NON_STATE ][ H11 ], i1->second[ NON_STATE ][ H11 ], decimal_nil, decimal_one);

		const decimal_t ibd_00 = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ IBD_STATE ][ H00 ], i1->second[ IBD_STATE ][ H00 ], decimal_nil, decimal_one);
		const decimal_t ibd_01 = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ IBD_STATE ][ H01 ], i1->second[ IBD_STATE ][ H01 ], decimal_nil, decimal_one);
		const decimal_t ibd_11 = decimal_err + approx<decimal_t>(k, i0->first, i1->first, i0->second[ IBD_STATE ][ H11 ], i1->second[ IBD_STATE ][ H11 ], decimal_nil, decimal_one);
		
		const decimal_t non_sum = non_00 + non_01 + non_11;
		const decimal_t ibd_sum = ibd_00 + ibd_01 + ibd_11;
		
		Model::emiss_type point;
		
		point[ NON_STATE ][ H00 ] = non_00 / non_sum;
		point[ NON_STATE ][ H01 ] = non_01 / non_sum;
		point[ NON_STATE ][ H11 ] = non_11 / non_sum;
		
		point[ IBD_STATE ][ H00 ] = ibd_00 / ibd_sum;
		point[ IBD_STATE ][ H01 ] = ibd_01 / ibd_sum;
		point[ IBD_STATE ][ H11 ] = ibd_11 / ibd_sum;
		
		this->input_emiss[ k ] = point;
	}
	
	
	// fill in for each site
	
	this->emiss.resize(grid->marker().size());
	
	Marker::Vector::const_iterator M, M_end = grid->marker().cend();
	
	for (M = grid->marker().cbegin(); M != M_end; ++M)
	{
		this->emiss[ M->index.value ] = this->input_emiss[ M->hap_count[H1] ];
	}
	
	// print
	if (file_ptr != nullptr)
	{
		this->print_emission(*file_ptr);
	}
	
	// finish
	this->emiss_good = true;
}


// expected emission probabilities

void LoadHMM::make_emission(const Gen::Grid::Data grid)
{
	static constexpr decimal_t allow = 0.0001;
	
	// get sample size of current dataset
	const size_t N = grid->sample().size() * 2; // number of sampled haplotypes
	
	const decimal_t n = static_cast<decimal_t>(N);
	
	for (size_t k = 1; k < N; ++k)
	{
		const decimal_t q = static_cast<decimal_t>(k) / n;
		const decimal_t p = decimal_one - q;
		
		Model::emiss_type point;
		
		point[ NON_STATE ][ H00 ] = std::pow(p, 2);
		point[ NON_STATE ][ H01 ] = 2 * p * q;
		point[ NON_STATE ][ H11 ] = std::pow(q, 2);
		
		decimal_t sum = p + q + (2 * p * q  *  allow);
		
		point[ IBD_STATE ][ H00 ] = p / sum;
		point[ IBD_STATE ][ H01 ] = (2 * p * q  *  allow) / sum;
		point[ IBD_STATE ][ H11 ] = q / sum;
		
		this->input_emiss[ k ] = point;
	}
	
	// fill in for each site
	
	this->emiss.resize(grid->marker().size());
	
	Marker::Vector::const_iterator M, M_end = grid->marker().cend();
	
	for (M = grid->marker().cbegin(); M != M_end; ++M)
	{
		this->emiss[ M->index.value ] = this->input_emiss[ M->hap_count[H1] ];
	}
	
	// finish
	this->emiss_good = true;
}


// generate transition probabilties

void LoadHMM::make_transition(const Gen::Grid::Data grid, std::ofstream * file_ptr)
{
	if (this->dists_good)
		throw std::runtime_error("Transition probabilities already generated");
	
	// get sample size of current dataset
	const size_t N = grid->sample().size() * 2; // number of sampled haplotypes
	
	if (N != this->Nh)
	{
		throw std::runtime_error("Unexpected sample size detected");
	}
	
	const size_t size = grid->marker().size();
	
	this->dists.resize(size - 1);
	
	for (size_t i = 1; i < size; ++i)
	{
		const decimal_t d0 = static_cast<decimal_t>(grid->marker(i-1).gen_dist);
		const decimal_t d1 = static_cast<decimal_t>(grid->marker( i ).gen_dist);
		
		const decimal_t dist = d1 - d0; // delta in cM
		
		// check negative
		if (dist < decimal_nil)
		{
			throw std::logic_error("Gen distance between positions is negative");
		}
		
		// insert into transition container
		
		this->dists[ i - 1 ] = dist + decimal_err;
	}
	
	// print
	if (file_ptr != nullptr)
	{
		this->print_transition(*file_ptr);
	}

	// finish
	this->dists_good = true;
}


// return model

Model::Data LoadHMM::make_model()
{
	if (this->done)
		throw std::runtime_error("Model already generated");
	
	if (!this->inits_good)
		throw std::runtime_error("Initial state probabilities not generated");
	if (!this->emiss_good)
		throw std::runtime_error("Emission probabilities not generated");
	if (!this->dists_good)
		throw std::runtime_error("Gen distances for transition probabilities not generated");
	
	this->done = true;
	
	return std::make_shared<Model>(this->Ne, this->Nh, std::move(this->inits_con), std::move(this->inits_dis), std::move(this->emiss), std::move(this->dists));
}


// print probabilities

void LoadHMM::print_initial(std::ofstream & stream) const
{
	if (!this->inits_good)
		throw std::runtime_error("Initial state probabilities not generated");
	
	// header
	stream << "Frequency" << ' ';
	stream << "CON_NON" << ' ';
	stream << "CON_IBD" << ' ';
	stream << "DIS_NON" << ' ';
	stream << "DIS_IBD" << std::endl;
	
	
	// content
	
	inits_input::const_iterator
	cit = this->input_inits_con.cbegin(),
	cti = this->input_inits_con.cend(),
	dit = this->input_inits_dis.cbegin(),
	dti = this->input_inits_dis.cend();
	
	while (cit != cti || dit != dti)
	{
		stream << std::fixed;
		stream << std::setprecision(8) << static_cast<decimal_t>(cit->first) / static_cast<decimal_t>(this->Nh) << ' ';
		stream << std::setprecision(8) << cit->second[NON_STATE] << ' ';
		stream << std::setprecision(8) << cit->second[IBD_STATE] << ' ';
		stream << std::setprecision(8) << dit->second[NON_STATE] << ' ';
		stream << std::setprecision(8) << dit->second[IBD_STATE] << std::endl;
		
		++cit;
		++dit;
	}
}

void LoadHMM::print_emission(std::ofstream & stream) const
{
	if (!this->emiss_good)
		throw std::runtime_error("Emission probabilities not generated");
	
	// header
	stream << "Frequency" << ' ';
	stream << "NON_00" << ' ';
	stream << "NON_01" << ' ';
	stream << "NON_11" << ' ';
	stream << "IBD_00" << ' ';
	stream << "IBD_01" << ' ';
	stream << "IBD_11" << std::endl;
	
	
	// content
	
	emiss_input::const_iterator it, ti = this->input_emiss.cend();
	
	for (it = this->input_emiss.cbegin(); it != ti; ++it)
	{
		stream << std::fixed;
		stream << std::setprecision(8) << static_cast<decimal_t>(it->first) / static_cast<decimal_t>(this->Nh) << ' ';
		
		stream << std::setprecision(8) << it->second[ NON_STATE ][ H00 ] << ' ';
		stream << std::setprecision(8) << it->second[ NON_STATE ][ H01 ] << ' ';
		stream << std::setprecision(8) << it->second[ NON_STATE ][ H11 ] << ' ';
		
		stream << std::setprecision(8) << it->second[ IBD_STATE ][ H00 ] << ' ';
		stream << std::setprecision(8) << it->second[ IBD_STATE ][ H01 ] << ' ';
		stream << std::setprecision(8) << it->second[ IBD_STATE ][ H11 ] << std::endl;
	}
}

void LoadHMM::print_transition(std::ofstream & stream) const
{
	if (!this->dists_good)
		throw std::runtime_error("Gen distances for transition probabilities not generated");
	
	// header
	stream << "Marker0" << ' ';
	stream << "Marker1" << ' ';
	stream << "Distance"  << std::endl;
	
	
	// content
	
	const size_t size = this->dists.size();
	
	for (size_t i = 0; i < size; ++i)
	{
		stream << i  << ' ';
		stream << i + 1 << ' ';
		stream << std::fixed << std::setprecision(8) << this->dists.at(i) << std::endl;
	}
}


void LoadHMM::print_initial(const std::string & filename) const
{
	std::ofstream stream(filename);
	this->print_initial(stream);
	stream.close();
}

void LoadHMM::print_emission(const std::string & filename) const
{
	std::ofstream stream(filename);
	this->print_emission(stream);
	stream.close();
}

void LoadHMM::print_transition(const std::string & filename) const
{
	std::ofstream stream(filename);
	this->print_transition(stream);
	stream.close();
}

