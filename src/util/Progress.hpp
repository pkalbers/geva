//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Progress_hpp
#define Progress_hpp


#include <math.h>
#include <stdio.h>

#include <chrono>
#include <mutex>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>



class Progress
{
public:
	
	// constructs
	Progress(const size_t, std::ostream & = std::cout); // progress bar, counting to target value
	Progress(std::ostream & = std::cout); // progress counter, unknown target value
	Progress(const size_t, const std::string &, std::ostream & = std::cout); // progress bar, with label
	Progress(const std::string &, std::ostream & = std::cout); // progress counter, with label
	
	// print new line
	void halt() const;
	
	// update iteration count
	void update(const size_t = 1);
	
	// set final iteration done
	void finish();
	
private:
	
	typedef std::chrono::steady_clock                        clock_t; // set clock type
	typedef std::chrono::milliseconds                         tick_t; // set clock resolution
	typedef std::chrono::time_point< clock_t, tick_t >        time_t; // set time measure type
	typedef std::chrono::duration< double, tick_t::period > period_t; // set duration type
	
	typedef std::chrono::duration< double, std::chrono::milliseconds::period > period_ms;
	typedef std::chrono::duration< double, std::chrono::seconds::period >      period_s;
	typedef std::chrono::duration< double, std::chrono::minutes::period >      period_m;
	typedef std::chrono::duration< double, std::chrono::hours::period >        period_h;
	
	typedef std::lock_guard<std::mutex> guard_t;
	
	
	// get current status
	std::string status(const int = 256) const;
	
	// print update message
	void print_update() const;
	
	// print final message
	void print_finish() const;
	
	const size_t iter; // expected end count
	const time_t time; // start time
	std::string  unit; // optional unit for labelling
	
	std::ostream & stream; // output stream
	
	time_t last; // time of last update
	size_t iter_i, // current count
		   skip_n, // number of skipped updates
		   skip_i; // countdown skipped updates
	
	std::mutex guard;
	
	// get difference between now and past timepoint
	static double interval(const time_t &);
};



#endif /* Progress_hpp */
