//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "LoadSim.hpp"


LoadSim::LoadSim(const std::string & filename)
: good(false)
{
	static const size_t header_fields = 8;
	static const std::string header[header_fields] = {
		"MarkerID",
		"SampleID0",
		"Chr0",
		"SampleID1",
		"Chr1",
		"Shared",
		"TrueLHS",
		"TrueRHS" };
	
	// open file
	this->stream.open(filename);
	
	// read first line
	if (this->stream.next())
	{
		Reader::Current line = this->stream.line();
		
		// count number of fields
		size_t n = 0;
		
		while (line.split())
		{
			Reader::Current field = line.field();
			
			if (field.str() != header[n])
			{
				throw std::invalid_argument("Unknown field: " + field.str());
			}
			
			++n;
		}
		
		if (n != header_fields)
		{
			throw std::invalid_argument("Unknown results format: " + filename);
		}
	}
	else
	{
		throw std::invalid_argument("Unable to read from results file: " + filename);
	}
	
	this->good = true;
}


// forward to next line

bool LoadSim::next()
{
	return (this->good && this->stream.next());
}


// parse current line

bool LoadSim::parse(size_t & site, IBD::SIM::Truth & element)
{
	Reader::Current line = this->stream.line();
	
	// walkabout
	while (line.split())
	{
		Reader::Current field = line.field();
		
		switch (field.number)
		{
			case 0: // MarkerID
			{
				site = field.convert<size_t>();
				break;
			}
			case 1: // SampleID0
			{
				element.pair.first = field.convert<size_t>();
				break;
			}
			case 2: // Chr0
			{
				const std::string chr = field.str();
				
				if (chr == "M") element.chr.first = Gen::ChrType::MATERNAL; else
				if (chr == "P") element.chr.first = Gen::ChrType::PATERNAL; else
					return false;
				
				break;
			}
			case 3: // SampleID1
			{
				element.pair.second = field.convert<size_t>();
				break;
			}
			case 4: // Chr1
			{
				const std::string chr = field.str();
				
				if (chr == "M") element.chr.second = Gen::ChrType::MATERNAL; else
				if (chr == "P") element.chr.second = Gen::ChrType::PATERNAL; else
					return false;
				
				break;
			}
			case 5: // Shared
			{
				int shr = field.convert<int>();
				
				if (shr == 0) element.shared = false; else
				if (shr == 1) element.shared = true; else
					return false;
				
				break;
			}
			case 6: // LHS
			{
				element.segment[IBD::LHS] = field.convert<size_t>();
				break;
			}
			case 7: // RHS
			{
				element.segment[IBD::RHS] = field.convert<size_t>();
				break;
			}
			default:
			{
				throw std::runtime_error("Unexpected element on line: " + std::to_string(line.number));
				break;
			}
		}
	}
	
	return true;
}

