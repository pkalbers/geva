//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef LoadGen_hpp
#define LoadGen_hpp


#include <stdio.h>
#include <exception>
#include <utility>
#include <set>
#include <sstream>

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenMap.hpp"
#include "GenGrid.hpp"

#include "Reader.hpp"


class LoadGen
{
public:
	
	// Filtering options
	struct Filter
	{
		int    chromosome      = -1;
		size_t position_beg    = 0;
		size_t position_end    = 0;
		bool   remove_missing  = true;
		bool   require_snp     = true;
	};
	
	
	// construct
	LoadGen(const std::string &, const std::string &);
	
	// forward to next line
	bool next();
	
	// parse current line
	bool parse(Gen::Grid::Make &, const Gen::Map &);
	
	// get detected chromosome (first occurrence is taken), -1 = unknown
	int chromosome() const;
	
	Gen::Sample::Vector sample; // sample vector
	Gen::Marker::Vector marker; // marker vector
	
	Filter filter; // filter settings
	
private:
	
	Reader stream; // stream of GEN file
	
	bool chrom_avail; // flag if chromsome was already scanned for
	int  chrom_value; // optional chromosome value
	
	bool good; // flag status, false = exit
};


#endif /* LoadGen_hpp */

