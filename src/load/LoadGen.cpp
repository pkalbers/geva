//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "LoadGen.hpp"


using namespace Gen;


LoadGen::LoadGen(const std::string & gen_file, const std::string & sample_file)
: chrom_avail(false)
, chrom_value(-1)
, good(false)
{
	// read sample file
	Reader sample_stream(sample_file);
	
	// skip header (2 lines)
	if (!sample_stream.next() ||
		!sample_stream.next())
	{
		throw std::runtime_error("Unable to read from SAMPLE file: " + sample_file);
	}
	
	std::set<std::string> unique;
	
	// walkabout samples
	while (sample_stream.next())
	{
		Reader::Current line = sample_stream.line();
		
		if (!line.split())
		{
			throw std::runtime_error("Invalid format of SAMPLE file: " + sample_file);
		}
		
		// parse sample
		Sample S;
		
		S.label = line.field().str();
		S.phase = false;
		
		unique.insert(S.label);
		this->sample.push_back(std::move(S));
		
		// check duplicates
		if (this->sample.size() != unique.size())
		{
			throw std::logic_error("Duplicate sample ID detected in SAMPLE file: " + sample_file);
		}
	}
	
	// open GEN file
	this->stream.open(gen_file);
	
	this->good = true;
}


// forward to next line

bool LoadGen::next()
{
	return (this->good && this->stream.next());
}


// parse current line

bool LoadGen::parse(Gen::Grid::Make & buffer, const Gen::Map & gmap)
{
	Reader::Current line = this->stream.line();
	
	int         chr = -1;
	std::string label0, label1;
	size_t      pos = 0, last_pos = 0;
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
			case 1: // LABEL 0
			{
				label0 = field.str();
				break;
			}
			case 2: // LABEL 1
			{
				label1 = field.str();
				break;
			}
			case 3: // POSITION
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
			case 4: // REF
			{
				ref = (field.size() == 1);
				allele += field.str();
				
				if (this->filter.remove_missing && field[0] == '.')
				{
					return false;
				}
				break;
			}
			case 5: // ALT
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
			default: // Genotypes
			{
				double g[3];
				
				unsigned i = 0;
				
				do
				{
					g[i++] = line.field().convert<double>();
					
					if (i == 3)
					{
						gen_t gt = index_to_genotype(G_);
						
						if (g[0] > g[1] && g[0] > g[2]) gt = index_to_genotype(G0); else
						if (g[1] > g[0] && g[1] > g[2]) gt = index_to_genotype(G1); else
						if (g[2] > g[0] && g[2] > g[1]) gt = index_to_genotype(G2);
						
						// insert genotype into buffer
						buffer.insert(gt);
						
						// allele and genotype counter
						M.count(gt);
						
						++count;

						i = 0;
					}
				}
				while (line.split());
				
				if (i != 0)
				{
					throw std::runtime_error("Incomplete number of genotypes on line " + std::to_string(line.number));
				}
				
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
	
	M.label = label0 + "_" + label1;
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

int LoadGen::chromosome() const
{
	if (!this->chrom_avail)
	{
		throw std::logic_error("Unable to return chromosome before reading from VCF file");
	}
	return this->chrom_value;
}

