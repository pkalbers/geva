//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#include "Progress.hpp"


#define PROGRESS_INTERVAL 2000 // milliseconds
#define PROGRESS_LINESIZE 80   // number of characters on printed line


// constructs

Progress::Progress(const size_t max_count, std::ostream & out_stream)
: iter(max_count)
, time(std::chrono::time_point_cast< tick_t >(clock_t::now()))
, stream(out_stream)
{
	this->last = this->time;
	this->iter_i = 0;
	this->skip_n = 1;
	this->skip_i = 1;
	
	this->stream << " ...\r" << std::flush;
}

Progress::Progress(std::ostream & out_stream)
: iter(0)
, time(std::chrono::time_point_cast< tick_t >(clock_t::now()))
, stream(out_stream)
{
	this->last = this->time;
	this->iter_i = 0;
	this->skip_n = 1;
	this->skip_i = 1;
	
	this->stream << " ...\r" << std::flush;
}

Progress::Progress(const size_t max_count, const std::string & label, std::ostream & out_stream)
: iter(max_count)
, time(std::chrono::time_point_cast< tick_t >(clock_t::now()))
, stream(out_stream)
{
	this->unit = (label.size() == 0) ? std::string(): label;
	
	this->last = this->time;
	this->iter_i = 0;
	this->skip_n = 1;
	this->skip_i = 1;
	
	this->stream << " ...\r" << std::flush;
}

Progress::Progress(const std::string & label, std::ostream & out_stream)
: iter(0)
, time(std::chrono::time_point_cast< tick_t >(clock_t::now()))
, stream(out_stream)
{
	this->unit = (label.size() == 0) ? std::string(): label;
	
	this->last = this->time;
	this->iter_i = 0;
	this->skip_n = 1;
	this->skip_i = 1;
	
	this->stream << " ...\r" << std::flush;
}


// print new line
void Progress::halt() const
{
	this->print_update();
	this->stream << std::endl;
}


// update iteration count
void Progress::update(const size_t step)
{
	guard_t lock(guard);
	
	this->iter_i += step;
	
	--this->skip_i;
	
	if (this->skip_i == 0)
	{
		this->skip_i = this->skip_n;
		
		if (Progress::interval(this->last) > PROGRESS_INTERVAL)
		{
			this->last = std::chrono::time_point_cast< tick_t >(clock_t::now());
			
			this->print_update();
		}
		else
		{
			++this->skip_n;
		}
	}
}


// set final iteration done
void Progress::finish()
{
	this->last = std::chrono::time_point_cast< tick_t >(clock_t::now());
	
	this->print_finish();
	
	this->stream << std::endl;
}


// get current status
std::string Progress::status(const int width) const
{
	static const double step_h = 1000 * 60 * 60;
	static const double step_m = 1000 * 60;
	static const double step_s = 1000;
	
	static const char * format_label_h  = " %lu %s in %.1f hours (%.1f / hour)";
	static const char * format_label_m  = " %lu %s in %.1f min (%.1f / min)";
	static const char * format_label_s  = " %lu %s in %.1f sec (%.1f / sec)";
	static const char * format_label_ms = " %lu %s in %.1f ms (%.1f / ms)";
	
	static const char * format_h  = " %lu in %.1f hours (%.1f / hour)";
	static const char * format_m  = " %lu in %.1f min (%.1f / min)";
	static const char * format_s  = " %lu in %.1f sec (%.1f / sec)";
	static const char * format_ms = " %lu in %.1f ms (%.1f / ms)";
	
	const bool no_label = this->unit.empty();
	
	const tick_t diff = (this->last - this->time);
	const char * form = (no_label) ? format_ms: format_label_ms;
	double       step = std::chrono::duration_cast< period_ms >(diff).count();
	
	char line[width];
	int  size = 0;
	
	// select timeframe
	if (step > step_h)
	{
		form = (no_label) ? format_h: format_label_h;
		step = std::chrono::duration_cast< period_h >(diff).count();
	}
	else if (step > step_m)
	{
		form = (no_label) ? format_m: format_label_m;
		step = std::chrono::duration_cast< period_m >(diff).count();
	}
	else if (step > step_s)
	{
		form = (no_label) ? format_s: format_label_s;
		step = std::chrono::duration_cast< period_s >(diff).count();
	}
	
	// calculate rate
	const double rate = (step <= 0.1) ? static_cast<double>(this->iter_i): static_cast<double>(this->iter_i) / step;
	
	// format with or w/out label
	if (no_label)
	{
		size = snprintf(line, width, form, this->iter_i, step, rate);
	}
	else
	{
		size = snprintf(line, width, form, this->iter_i, this->unit.c_str(), step, rate);
	}
	
	// format to string
	
	if (size < 1)
	{
		return " ... " + std::string(width - 5, ' ');
	}
	
	if (size > width)
	{
		line[width - 3] = '\0';
		return std::string(line) + "...";
	}
	
	return std::string(line) + std::string(width - size, ' ');
}


// print update message
void Progress::print_update() const
{
	static constexpr double width = static_cast<double>(PROGRESS_LINESIZE - 10);
	
	if (this->iter < this->iter_i)
	{
		this->stream << this->status(PROGRESS_LINESIZE - 1) << '\r' << std::flush;
	}
	else
	{
		const double ratio = static_cast<double>(this->iter_i) / static_cast<double>(this->iter);
		const size_t width_done = static_cast<size_t>(width * ratio + 0.5);
		const size_t width_rest = static_cast<size_t>(width - width_done);
		
		// print to stream
		this->stream
		<< std::setw(5) << std::right << std::setprecision(1) << std::fixed << (ratio * 100.0)
		<< "% [" << std::string(width_done, '=') << std::string(width_rest, ' ') << "]\r"
		<< std::flush;
	}
}


// print final message
void Progress::print_finish() const
{
	this->stream << "Completed" << this->status() << std::endl;
}


// get difference between now and past timepoint
double Progress::interval(const time_t & timepoint)
{
	const time_t now = std::chrono::time_point_cast< tick_t >(clock_t::now());
	
	return period_t(now - timepoint).count();
}

