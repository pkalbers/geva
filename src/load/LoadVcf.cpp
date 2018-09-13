//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//


#include "LoadVcf.hpp"


using namespace Gen;


LoadVcf::LoadVcf(const std::string & filename)
: chrom_avail(false)
, chrom_value(-1)
, good(false)
{
	static const size_t header_fields = 9;
	static const std::string header[header_fields] = {
		"#CHROM",
		"POS",
		"ID",
		"REF",
		"ALT",
		"QUAL",
		"FILTER",
		"INFO",
		"FORMAT" };
	
	this->stream.open(filename);
	
	// read first line
	if (this->stream.next())
	{
		Reader::Current line = this->stream.line();
		
		// check format on first line
		if (line.size() < 16 || std::string(line, 16) != "##fileformat=VCF")
		{
			throw std::logic_error("Input file not in Variant Call Format: " + filename);
		}
	}
	else
	{
		throw std::logic_error("Unable to read from VCF file: " + filename);
	}
	
	
	// walkabout header
	while(this->stream.next())
	{
		Reader::Current line = this->stream.line();
		
		// check if last header line reached
		if (line.size() > 6 && std::string(line, 6) == header[0])
		{
			std::ostringstream oss;
			bool   error = false;
			bool   check[ header_fields ];
			size_t i = 0;
			
			// look for required fields
			
			while (line.split())
			{
				check[i] = (header[i] == line.field().str());
				
				if (++i == header_fields)
				{
					break;
				}
			}
			
			for (i = 0; i < header_fields; ++i)
			{
				if (!check[i])
				{
					oss << std::endl << ' ' << header[i];
					error = true;
				}
			}
			
			if (error)
			{
				throw std::logic_error("Missing fields in VCF file: " + filename + oss.str());
			}
			
			// parse sample
			
			std::set<std::string> unique;
			
			while (line.split())
			{
				Sample S;
				
				S.label = line.field().str();
				S.phase = true;
				
				unique.insert(S.label);
				this->sample.push_back(std::move(S));
			}
			
			if (this->sample.size() != unique.size())
			{
				throw std::logic_error("Duplicate sample ID detected in VCF file: " + filename);
			}
			
			break;
		}
		
		// check if end of header reached before stopping
		if (line[0] != '#')
		{
			throw std::logic_error("Invalid header detected in VCF file: " + filename);
		}
	}
	
	this->good = true;
}


// forward to next line

bool LoadVcf::next()
{
	return (this->good && this->stream.next());
}


// parse current line

bool LoadVcf::parse(Gen::Grid::Make & buffer, const Gen::Map & gmap)
{
	Reader::Current line = this->stream.line();
	
	int         chr = -1;
	size_t      pos = 0, last_pos = 0;
	std::string label;
	bool        ref = false, alt = false;
	std::string allele;
	size_t count = 0;
	
	allele.reserve(256);
	
	Marker M;
	
	// walkabout
	while (line.split())
	{
		Reader::Current field = line.field();
		
		switch (field.number)
		{
			case 0: // CHROM
			{
				chr = field.convert<int>();
				
				if (this->filter.chromosome != -1 && this->filter.chromosome != chr)
				{
					return false;
				}
				
				if (this->chrom_avail)
				{
					if (this->chrom_value != chr)
					{
						return false;
					}
				}
				else
				{
					this->chrom_avail = true;
					this->chrom_value = chr;
				}
				break;
			}
			case 1: // POSITION
			{
				pos = field.convert<size_t>();
				
				if (this->filter.position_beg < this->filter.position_end)
				{
					if (pos < this->filter.position_beg)
					{
						return false;
					}
					if (pos >= this->filter.position_end)
					{
						this->good = false;
						return false;
					}
				}
				
				if (pos <= last_pos)
				{
					throw std::string("Invalid position on line " + std::to_string(line.number));
				}
				
				last_pos = pos;
				break;
			}
			case 2: // ID
			{
				label = field.str();
				break;
			}
			case 3: // REF
			{
				ref = (field.size() == 1);
				allele += field.str();
				
				if (this->filter.remove_missing && field[0] == '.')
				{
					return false;
				}
				break;
			}
			case 4: // ALT
			{
				alt = (field.size() == 1);
				allele += ',';
				allele += field.str();
				
				if (this->filter.remove_missing && field[0] == '.')
				{
					return false;
				}
				if (this->filter.require_snp && !(ref && alt))
				{
					return false;
				}
				break;
			}
			case 5: // QUAL
			{
				if (this->filter.require_qual_above > 0 && this->filter.require_qual_above > field.convert<int>())
				{
					return false;
				}
				break;
			}
			case 6: // FILTER
			{
				if (this->filter.require_filter_pass && field.str() != "PASS") // check filtering passed
				{
					return false;
				}
				break;
			}
			case 7: // INFO
			{
				// ignore
				break;
			}
			case 8: // FORMAT
			{
				bool flag = true;
				for (int i = 0, end = int(field.size()) - 1; i < end; ++i)
				{
					if (field[i] == 'G' && field[i+1] == 'T') // search "GT"
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					throw std::string("Missing GT format on line " + std::to_string(line.number));
				}
				break;
			}
			default: // Genotypes
			{
				do
				{
					Reader::Current g = line.field();
					
					if (g.size() < 3 || (g[1] != '/' && g[1] != '|'))
					{
						throw std::runtime_error("Invalid genotype on line " + std::to_string(line.number) + ", field " + std::to_string(g.number) + ":\n" + g.str());
					}
					
					const char c0 = g[0];
					const char c1 = g[2];
					const bool ph = (g[1] == '|');
					
					const gen_t gt = make_genotype(c0, c1, ph);
					
					// allele and genotype counter
					M.count(gt);
					
					// insert genotype into buffer
					buffer.insert(gt);
					
					// determine phasing
					if (!is_genotype<G_>(gt) && !ph && this->sample.at(count).phase)
					{
						this->sample.at(count).phase = false;
					}
					
					
					++count;
				}
				while (line.split());
				break;
			}
		}
	}
	
	if (count != this->sample.size())
	{
		throw std::runtime_error("Unexpected number of genotypes on line " + std::to_string(line.number) +
								 "\n Expected: " + std::to_string(this->sample.size()) +
								 "\n Detected: " + std::to_string(count));
	}
	
	if (!buffer.good)
	{
		throw std::runtime_error("Unexpected buffer error");
	}
	
	
	M.label = std::move(label);
	M.chromosome = chr;
	M.position   = pos;
	M.allele.parse(allele);
	
	
	// Approximate rate and distance
	
	Map::Element mapped = gmap.get(chr, pos);
	
	if (!mapped.valid())
	{
		throw std::invalid_argument("Invalid genetic map");
	}
	
	M.rec_rate = mapped.rate;
	M.gen_dist = mapped.dist;
	
	this->marker.push_back(std::move(M));
	
	return true;
}


// get detected chromosome (first occurrence is taken), -1 = unknown

int LoadVcf::chromosome() const
{
	if (!this->chrom_avail)
	{
		throw std::logic_error("Unable to return chromosome before reading from VCF file");
	}
	return this->chrom_value;
}

