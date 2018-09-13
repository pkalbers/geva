//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//


#include "GenGrid.hpp"


using namespace Gen;


// constant char sequence as identifier for binary file navigation
const char Grid::checkpoint[4] = { 0x25, 0x50, 0x4b, 0x41 };


// constructs

Grid::Grid(const std::string & filename)
: source(filename, Binary::mode::READ)
{
	this->load();
	this->load_sample();
	this->load_marker();
	
	this->buffer_limit = 0;
	this->buffer_count = 0;
}

Grid::Grid(Grid && other) // move
: size_sample(other.size_sample)
, size_marker(other.size_marker)
, interval(std::move(other.interval))
, compression(other.compression)
, source(std::move(other.source))
, buffer(std::move(other.buffer))
, buffer_limit(other.buffer_limit)
, buffer_count(other.buffer_count)
, sample_list(std::move(other.sample_list))
, marker_list(std::move(other.marker_list))
{}

Grid::Grid(Binary && bin)
: source(std::move(bin))
{
	this->load();
}


// read from file

Grid::Vector Grid::read(const Sample::Key & key, const bool decompress)
{
	guard_t lock(this->guard);
	
	this->source.jump(key);
	this->source.match<char>(checkpoint, 4);
	this->source.match<Bin4>(key.value);
	
	const size_t out_length = this->source.read<Bin4, size_t>(); // marker size
	const size_t raw_length = this->source.read<Bin4, size_t>(); // length
	
	Vector raw(raw_length);
	
	this->source.read<Bin1>(raw.data(), raw_length); // data vector
	
	if (this->compression && decompress)
	{
		Vector out = decompress_genotype_vector(raw, raw_length, out_length);
		
		if (out.size() != out_length)
			throw std::runtime_error("Unable to decompress data vector");
		
		return out;
	}
	
	return raw;
}


// fetch from cache

Variant::Vector::Data Grid::get(const Sample::Key & key)
{
	guard_t lock(this->guard);
	
	Variant::Vector::Data & ptr = this->buffer.at(key.value);
	
	if (ptr)
	{
		return ptr;
	}
	
	this->prune();
	
	this->source.jump(key);
	this->source.match<char>(checkpoint, 4);
	this->source.match<Bin4>(key.value);
	
	const size_t out_length = this->source.read<Bin4, size_t>(); // marker size
	const size_t raw_length = this->source.read<Bin4, size_t>(); // length
	
	Vector raw(raw_length);
	
	this->source.read<Bin1>(raw.data(), raw_length); // data vector
	
	if (this->compression)
	{
		Vector out = decompress_genotype_vector(raw, raw_length, out_length);
		
		if (out.size() != out_length)
			throw std::runtime_error("Unable to decompress data vector");
		
		ptr = std::make_shared< Variant::Vector >(std::move(out), this->sample_list[ key.value ].phase);
	}
	else
	{
		ptr = std::make_shared< Variant::Vector >(std::move(raw), this->sample_list[ key.value ].phase);
	}

	++this->buffer_count;
	
	return ptr;
}


// limit cache size

void Grid::cache(const size_t max)
{
	guard_t lock(this->guard);
	this->buffer_limit = max;
	this->prune();
}


// prune buffer

void Grid::prune()
{
	if (this->buffer_count > this->buffer_limit)
	{
		const size_t end = random_number(this->size_sample);
		
		size_t cyc = end + 1;
		
		while (cyc != end)
		{
			if (cyc == this->size_sample)
				cyc = 0;
			else
				++cyc;
			
			if (this->buffer[cyc].unique())
			{
				this->buffer[cyc].reset();
				--this->buffer_count;
			}
			
			if (this->buffer_count == this->buffer_limit)
				break;
		}
	}
}


// get sample/marker information

Sample::Reference Grid::sample(const Sample::Key & index) const
{
	return this->sample_list.at(index);
}

Marker::Reference Grid::marker(const Marker::Key & index) const
{
	return this->marker_list.at(index);
}


// get direct reference to sample/marker vector

Sample::Vector const & Grid::sample() const
{
	return this->sample_list;
}

Marker::Vector const & Grid::marker() const
{
	return this->marker_list;
}


// print samples/markers

void Grid::print_sample(std::ostream & stream) const
{
	// header
	stream << Sample::header << std::endl;
	
	// content
	
	Sample::Vector::const_iterator it, ti = this->sample_list.cend();
	
	for (it = this->sample_list.cbegin(); it != ti; ++it)
	{
		it->print(stream);
		stream << std::endl;
	}
}

void Grid::print_marker(std::ostream & stream) const
{
	// header
	stream << Marker::header << std::endl;
	
	// content
	
	Marker::Vector::const_iterator it, ti = this->marker_list.cend();
	
	for (it = this->marker_list.cbegin(); it != ti; ++it)
	{
		it->print(stream);
		stream << std::endl;
	}
}

void Grid::print_sample(const std::string & filename) const
{
	std::ofstream stream(filename, std::fstream::out);
	this->print_sample(stream);
	stream.close();
}

void Grid::print_marker(const std::string & filename) const
{
	std::ofstream stream(filename, std::fstream::out);
	this->print_marker(stream);
	stream.close();
}


// read from file

void Grid::load()
{
	// header
	this->source.match<char>(checkpoint, 4);
	this->source.read<Bin4>(this->size_sample); // sample size
	this->source.read<Bin4>(this->size_marker); // marker size
	this->source.read<Bin4>(this->interval.data(), this->interval.size()); // current marker interval
	this->source.read<bool>(this->compression); // compression setting
	
	// walkabout content
	for (size_t i = 0; i < this->size_sample; ++i)
	{
		// mark position
		this->source.here(i);
		
		// check index and marker size
		this->source.match<char>(checkpoint, 4);
		const size_t key = this->source.read<Bin4, size_t>();
		const size_t out_length = this->source.read<Bin4, size_t>();
		const size_t raw_length = this->source.read<Bin4, size_t>();
		
		if (key != i || out_length != this->size_marker)
		{
			throw std::invalid_argument("Invalid binary file format");
		}
		
		// skip length of vector
		this->source.skip<Type>(raw_length);
		
		// set buffer
		this->buffer[ i ] = nullptr;
	}
	
	// separator
	this->source.match<char>(checkpoint, 4);
}


// read sample/marker information

void Grid::load_sample()
{
	this->sample_list.reserve(this->size_sample);
	
	for (size_t i = 0; i < this->size_sample; ++i)
	{
		Sample S;
		size_t n;
		
		this->source.match<char>(checkpoint, 4);
		this->source.read<Bin4>(S.index.value); // index
		this->source.read<Bin4>(n); // length of label
		std::vector<char> label(n, '\0');
		this->source.read<char>(label.data(), n); // label
		S.label = std::string(label.data(), label.size());
		this->source.read<bool>(S.phase); // phase
		
		if (S.index.value != i)
		{
			throw std::runtime_error("Unexpected error while loading sample information");
		}
		
		this->sample_list.push_back(std::move(S));
	}
	
	this->source.match<char>(checkpoint, 4);
}

void Grid::load_marker()
{
	this->marker_list.reserve(this->size_marker);
	
	for (size_t i = 0; i < this->size_marker; ++i)
	{
		Marker M;
		size_t n;
		
		this->source.match<char>(checkpoint, 4);
		this->source.read<Bin4, size_t>(M.index.value); // index
		this->source.read<Bin4>(n); // length of label
		std::vector<char> label(n, '\0');
		this->source.read<char>(label.data(), n); // label
		M.label = std::string(label.data(), label.size());
		this->source.read<Bin2>(M.chromosome); // chromosome
		this->source.read<Bin4>(M.position); // position
		this->source.read<Bin4>(n); // length of allele
		std::vector<char> allele(n, '\0');
		this->source.read<char>(allele.data(), n); // allele
		M.allele.parse(std::string(allele.data(), allele.size()));
		this->source.read<Bin4>(M.hap_count.data(), M.hap_count.size()); // allele count
		this->source.read<Bin4>(M.gen_count.data(), M.gen_count.size()); // genotype count
		this->source.read<double>(M.rec_rate); // recombination rate
		this->source.read<double>(M.gen_dist); // genetic distance
		
		if (M.index.value != i)
		{
			throw std::runtime_error("Unexpected error while loading marker information");
		}
		
		this->marker_list.push_back(std::move(M));
	}
	
	this->source.match<char>(checkpoint, 4);
}



//
// Make new grid
//

Grid::Make::Make(const std::string & filename, const size_t & sample_size, const size_t & marker_size, const bool compress_data)
: unique(random_string(16))
, output(filename)
, comprs(compress_data)
, matrix(matrix_t(sample_size, vector_t(marker_size)))
{
	this->irow = 0;
	this->icol = 0;
	this->nrow = sample_size;
	this->ncol = marker_size;
	
	this->interval[0] = 0;
	this->interval[1] = 0;
	
	this->good = false;
	this->full = false;
}


// fill buffer with genotypes

void Grid::Make::insert(const gen_t g)
{
	if (this->full)
	{
		throw std::runtime_error("Cannot inset into full Buffer");
	}
	
	// insert genotype in buffer matrix
	this->matrix[ this->irow++ ][ this->icol ] = g;
	
	this->good = false;
	
	// sample complete
	if (this->irow == this->nrow)
	{
		this->irow = 0;
		
		// count marker
		++this->icol;
		++this->interval[1];
		
		this->good = true;
		
		// matrix filled
		if (this->icol == this->ncol)
		{
			this->full = true;
		}
	}
}


// write buffer to binary file and reset buffer

void Grid::Make::save(const bool auto_delete)
{
	if (!this->good)
	{
		throw std::runtime_error("Unwilling to save incomplete marker vector");
	}
	
	// create temporary file
	if (auto_delete)
	{
		std::ostringstream tmp;
		
		tmp << "tmp" << '.';
		tmp << std::setw(4) << std::setfill('0') << this->sources.size() << '.';
		tmp << this->unique << '.';
		tmp << random_string(16) << '.';
		tmp << "bin";
		
		const std::string tmp_filename = tmp.str();
		
		this->sources.push_back(Binary(tmp_filename, Binary::mode::WRITE, true));
	}
	else
	{
		this->sources.push_back(Binary());
	}
	
	source_t::reference bin = this->sources.back();
	
	
	// header
	bin.write<char>(checkpoint, 4);
	bin.write<Bin4>(this->nrow); // sample size
	bin.write<Bin4>(this->icol); // marker size
	bin.write<Bin4>(this->interval.data(), this->interval.size()); // current marker interval
	bin.write<bool>(this->comprs); // compression setting
	
	
	if (this->comprs) // compression
	{
		// walkabout samples
		for (this->irow = 0; this->irow < this->nrow; ++this->irow)
		{
			const vector_t cmprss = compress_genotype_vector(this->matrix[ this->irow ], this->icol);
			const size_t   length = cmprss.size();

			// write vector
			bin.write<char>(checkpoint, 4);
			bin.write<Bin4>(this->irow); // index
			bin.write<Bin4>(this->icol); // marker size
			bin.write<Bin4>(length); // length
			bin.write<Bin1>(cmprss.data(), length); // data vector
			
			// reset vector
			this->matrix[ this->irow ] = vector_t(this->ncol);
		}
	}
	else
	{
		// walkabout samples
		for (this->irow = 0; this->irow < this->nrow; ++this->irow)
		{
			// write vector
			bin.write<char>(checkpoint, 4);
			bin.write<Bin4>(this->irow); // index
			bin.write<Bin4>(this->icol); // marker size
			bin.write<Bin4>(this->icol); // length
			bin.write<Bin1>(this->matrix[ this->irow ].data(), this->icol); // data vector
			
			// reset vector
			this->matrix[ this->irow ] = vector_t(this->ncol);
		}
	}
	
	
	// footer
	bin.write<char>(checkpoint, 4); // EOF
	
	
	// reset buffer
	
	this->irow = 0;
	this->icol = 0;
	
	this->interval[0] = this->interval[1];
	
	this->good = false;
	this->full = false;
}


// save overhang and concatenate source files

void Grid::Make::finish(Sample::Vector & samples, Marker::Vector & markers)
{
	// save overhang
	if (this->good)
	{
		this->save();
	}
	else if (this->irow > 0)
	{
		throw std::runtime_error("Incomplete marker vector");
	}
	
	
	// load temporary files into grid
	std::vector<Grid> grids;
	
	interval_t full_interval{0, 0};
	
	source_t::iterator it, end = this->sources.end();
	
	for (it = this->sources.begin(); it != end; ++it)
	{
		it->reset(); // rewind file
		
		grids.push_back(Grid(std::move(*it))); // move binary file
		
		std::vector<Grid>::const_reference grid = grids.back();
		
		// check sample size
		if (samples.size() != grid.size_sample)
		{
			throw std::runtime_error("Unexpected error while reading temporary file");
		}
		
		// get full marker interval
		if (full_interval[0] > grid.interval[0]) full_interval[0] = grid.interval[0];
		if (full_interval[1] < grid.interval[1]) full_interval[1] = grid.interval[1];
	}
	this->sources.clear(); // clear binary sources
	
	
	// marker size
	size_t full_length = full_interval[1] - full_interval[0];
	
	
	// check sample and marker size
	if (samples.size() != this->nrow || markers.size() != full_length)
	{
		throw std::runtime_error("Unexpected error while reading temporary file");
	}
	
	
	// create output file
	Binary out(this->output, Binary::mode::WRITE);
	
	
	// write header
	out.write<char>(checkpoint, 4);
	out.write<Bin4>(this->nrow); // sample size
	out.write<Bin4>(full_length); // marker size
	out.write<Bin4>(full_interval.data(), full_interval.size()); // current marker interval
	out.write<bool>(this->comprs); // compression setting
	
	
	// walkabout samples
	std::vector<Grid>::iterator grid, last = grids.end();
	
	for (size_t i = 0; i < this->nrow; ++i)
	{
		const Sample::Key key(i);
		
		vector_t full;
		full.reserve(full_length);
		
		// walkabout grids
		for (grid = grids.begin(); grid != last; ++grid)
		{
			Vector part = grid->read(key, false); // fetch part length vector
			
			full.insert(full.end(), part.begin(), part.end());
		}
		
		const size_t size = full.size();
		
		if (!this->comprs && size != full_length)
		{
			throw std::runtime_error("Temporary file was corrupted");
		}
		
		// write data vector
		out.write<char>(checkpoint, 4);
		out.write<Bin4>(i); // index
		out.write<Bin4>(full_length); // marker size
		out.write<Bin4>(size); // length
		out.write<Bin1>(full.data(), size); // data vector
	}
	
	
	// write separator
	out.write<char>(checkpoint, 4);
	
	
	// write sample information
	for (size_t i = 0; i < this->nrow; ++i)
	{
		samples[i].index = i; // assign index
		save_sample(out, samples[i]);
	}
	
	// write separator
	out.write<char>(checkpoint, 4);
	
	
	// write marker information
	for (size_t i = 0; i < full_length; ++i)
	{
		markers[i].index = i; // assign index
		save_marker(out, markers[i]);
	}
	
	// write footer
	out.write<char>(checkpoint, 4);
}



// Join multiple grids

// Grid source

// construct

Grid::Join::Source::Source(const Grid::Data ptr)
: grid(ptr)
, position_min(ptr->marker().front().position)
, position_max(ptr->marker().back().position)
{}


// sort, check overlap

bool Grid::Join::Source::operator <  (const Grid::Join::Source & other) const
{
	return (this->position_min < other.position_min);
}
bool Grid::Join::Source::operator  > (const Grid::Join::Source & other) const
{
	return (this->position_min > other.position_min);
}
bool Grid::Join::Source::operator == (const Grid::Join::Source & other) const
{
	return ((other.position_min <= this->position_min && this->position_min <= other.position_max ) ||
			(other.position_min <= this->position_max && this->position_max <= other.position_max ));
}


// construct

Grid::Join::Join(const std::string & filename)
: output(filename)
, compression(false)
, marker_size(0)
, sample_size(0)
, good(false)
{}


// insert grid

void Grid::Join::insert(const Grid::Data input)
{
	Source join(input);
	
	if (this->good)
	{
		if (this->compression != join.grid->compression)
			throw std::invalid_argument("Input files have different compression setting");
		
		if (this->sample_size != join.grid->sample_size())
			throw std::invalid_argument("Input files have different sample sizes");
		
		
		// compare samples
		
		for (size_t i = 0; i < this->sample_size; ++i)
		{
			if (this->sample_list[i] != join.grid->sample(i))
				throw std::invalid_argument("Input files have different samples");
		}
		
		
		// check overlap
		
		source_t::const_iterator it, ti = this->source.cend();
		
		for (it = this->source.cbegin(); it != ti; ++it)
		{
			if (*it == join)
				throw std::invalid_argument("Input files have overlapping markers");
		}
	}
	else
	{
		this->compression = join.grid->compression;
		this->sample_size = join.grid->sample_size();
		this->sample_list = join.grid->sample_list;
		
		this->good = true;
	}
	
	// count markers
	this->marker_size += join.grid->marker_size();
	
	// append grid
	this->source.push_back(join);
}


// join inserted grids

void Grid::Join::finish()
{
	if (!this->good)
	{
		throw std::runtime_error("No input files provided");
	}
	
	
	// sort grid list
	this->source.sort();
	
	source_t::const_iterator join, join_end = this->source.cend();
	
	
	// join marker lists
	
	Marker::Vector marker_list;
	Grid::interval_t  interval;
	
	interval[0] = 0;
	interval[1] = 0;
	
	for (join = this->source.cbegin(); join != join_end; ++join)
	{
		interval[1] += join->grid->marker_size();
		
		marker_list.insert(marker_list.end(), join->grid->marker_list.begin(), join->grid->marker_list.end());
	}
	
	
	// check joined marker list
	
	size_t last_pos = 0;
	
	Marker::Vector::const_iterator marker, marker_end = marker_list.cend();
	
	for (marker = marker_list.cbegin(); marker != marker_end; ++marker)
	{
		if (marker->position < last_pos)
			throw std::runtime_error("Invalid marker vector in joined input files");
		
		last_pos = marker->position;
	}
	
	
	// create output file
	Binary out(this->output, Binary::mode::WRITE);
	
	
	// write header
	out.write<char>(checkpoint, 4);
	out.write<Bin4>(this->sample_size); // sample size
	out.write<Bin4>(this->marker_size); // marker size
	out.write<Bin4>(interval.data(), interval.size()); // current marker interval
	out.write<bool>(this->compression); // compression setting
	
	
	// walkabout samples
	
	for (size_t i = 0; i < this->sample_size; ++i)
	{
		const Sample::Key key(i);
		
		vector_t full;
		full.reserve(this->marker_size);
		
		// walkabout grids
		for (join = this->source.cbegin(); join != join_end; ++join)
		{
			Vector part = join->grid->read(key, false); // fetch part length vector
			
			full.insert(full.end(), part.begin(), part.end());
		}
		
		const size_t size = full.size();
		
		if (!this->compression && size != this->marker_size)
		{
			throw std::runtime_error("Input file was corrupted");
		}
		
		// write data vector
		out.write<char>(checkpoint, 4);
		out.write<Bin4>(i); // index
		out.write<Bin4>(this->marker_size); // marker size
		out.write<Bin4>(size); // length
		out.write<Bin1>(full.data(), size); // data vector
	}
	
	// write separator
	out.write<char>(checkpoint, 4);
	
	
	// write sample information
	for (size_t i = 0; i < this->sample_size; ++i)
	{
		this->sample_list[i].index = i; // assign index
		save_sample(out, this->sample_list[i]);
	}
	
	// write separator
	out.write<char>(checkpoint, 4);
	
	
	// write marker information
	for (size_t i = 0; i < this->marker_size; ++i)
	{
		marker_list[i].index = i; // assign index
		save_marker(out, marker_list[i]);
	}
	
	
	// write footer
	out.write<char>(checkpoint, 4);
}



// save sample/marker information

void Grid::save_sample(Binary & bin, const Sample & sample)
{
	bin.write<char>(checkpoint, 4);
	bin.write<Bin4>(sample.index.value); // index
	bin.write<Bin4>(sample.label.size()); // length of label
	bin.write<char>(sample.label.data(), sample.label.size()); // label
	bin.write<bool>(sample.phase); // phase
}

void Grid::save_marker(Binary & bin, const Marker & marker)
{
	const std::string allele = marker.allele.str();
	
	bin.write<char>(checkpoint, 4);
	bin.write<Bin4>(marker.index.value); // index
	bin.write<Bin4>(marker.label.size()); // length of label
	bin.write<char>(marker.label.data(), marker.label.size()); // label
	bin.write<Bin2>(marker.chromosome); // chromosome
	bin.write<Bin4>(marker.position); // position
	bin.write<Bin4>(allele.size()); // length of allele
	bin.write<char>(allele.data(), allele.size()); // allele
	bin.write<Bin4>(marker.hap_count.data(), marker.hap_count.size()); // allele count
	bin.write<Bin4>(marker.gen_count.data(), marker.gen_count.size()); // genotype count
	bin.write<double>(marker.rec_rate); // recombination rate
	bin.write<double>(marker.gen_dist); // genetic distance
}

