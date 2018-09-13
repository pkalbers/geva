//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef LoadVcf_hpp
#define LoadVcf_hpp


#include <exception>
#include <deque>
#include <utility>
#include <set>
#include <sstream>

#include "Gen.hpp"
#include "GenSample.hpp"
#include "GenMarker.hpp"
#include "GenMap.hpp"
#include "GenGrid.hpp"

#include "Reader.hpp"


class LoadVcf
{
public:
	
	// Filtering options
	struct Filter
	{
		int    chromosome            = -1;
		size_t position_beg          = 0;
		size_t position_end          = 0;
		bool   require_snp           = true;
		bool   remove_missing        = true;
		int    require_qual_above    = -1;
		bool   require_filter_pass   = false;
	};
	
	
	// construct
	LoadVcf(const std::string &);
	
	// forward to next line
	bool next();
	
	// parse current line
	bool parse(Gen::Grid::Make &, const Gen::Map &);
	
	// get detected chromosome (first occurrence is taken), -1 = unknown
	int chromosome() const;
	
	
	Gen::Marker::Vector marker; // marker vector
	Gen::Sample::Vector sample; // sample vector
	
	Filter filter; // filter settings
	

private:
	
	Reader stream; // stream of VCF file
	
	bool chrom_avail; // flag if chromsome was already scanned for
	int  chrom_value; // optional chromosome value
	
	bool good; // flag status, false = exit
};


#endif /* LoadVcf_hpp */

