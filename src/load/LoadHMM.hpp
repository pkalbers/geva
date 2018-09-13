//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef LoadHMM_hpp
#define LoadHMM_hpp

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "Approx.h"
#include "Decimal.h"

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenGrid.hpp"
#include "GenShare.hpp"

#include "IBD_HMM.hpp"

#include "Reader.hpp"


// Load HMM input data
class LoadHMM
{
public:
	
	// construct
	LoadHMM(const size_t &, const size_t &);
	
	// read initial state probabilties from file and generate matrix
	void make_initial(const Gen::Grid::Data, const std::string &, std::ofstream * = nullptr);
	
	// expected initial probabilities
	void make_initial(const Gen::Grid::Data);
	
	// read emission probabilties from file and generate matrix
	void make_emission(const Gen::Grid::Data, const std::string &, std::ofstream * = nullptr);
	
	// expected emission probabilities
	void make_emission(const Gen::Grid::Data);
	
	// prepare transition probabilties
	void make_transition(const Gen::Grid::Data, std::ofstream * = nullptr);
	
	// return model
	IBD::HMM::Model::Data make_model();
	
	// print probabilities
	void print_initial(std::ofstream &) const;
	void print_emission(std::ofstream &) const;
	void print_transition(std::ofstream &) const;
	void print_initial(const std::string &) const;
	void print_emission(const std::string &) const;
	void print_transition(const std::string &) const;
	
	
private:
	
	using inits_input = std::map< size_t, IBD::HMM::Model::inits_type >;
	using emiss_input = std::map< size_t, IBD::HMM::Model::emiss_type >;
	
	
	const size_t Ne;
	const size_t Nh;
	
	IBD::HMM::Model::inits_list inits_con; // concordant initial probabilties
	IBD::HMM::Model::inits_list inits_dis; // discordant initial probabilties
	IBD::HMM::Model::emiss_list emiss; // emission probabilties
	IBD::HMM::Model::dists_list dists; // distances, to prepare transition probabilties
	
	
	// internal
	
	inits_input input_inits_con;
	inits_input input_inits_dis;
	emiss_input input_emiss;
	
	bool inits_good;
	bool emiss_good;
	bool dists_good;
	bool done;	
};


#endif /* LoadHMM_hpp */

