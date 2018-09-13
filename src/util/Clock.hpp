//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Clock_hpp
#define Clock_hpp


// C headers (sorted)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// C++ headers
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

/*
namespace Time
{

	using clock_measure_t = double;
	
	using clock_process_clock = std::chrono::system_clock;
	using clock_elasped_clock = std::chrono::steady_clock;
	
	using clock_unit_n = std::chrono::nanoseconds;
	using clock_unit_u = std::chrono::milliseconds;
	using clock_unit_s = std::chrono::seconds;
	using clock_unit_m = std::chrono::minutes;
	using clock_unit_h = std::chrono::hours;
	
	using clock_duration_u = std::chrono::duration< clock_measure_t, clock_unit_u::period >;
	using clock_duration_s = std::chrono::duration< clock_measure_t, clock_unit_s::period >;
	using clock_duration_m = std::chrono::duration< clock_measure_t, clock_unit_m::period >;
	using clock_duration_h = std::chrono::duration< clock_measure_t, clock_unit_h::period >;
	
	using clock_resolution_t = clock_unit_u;

	
	template< class T >
	using clock_point_t = std::chrono::time_point< T, clock_resolution_t >;
	
	
	// get current time point
	template< class T >
	constexpr clock_point_t<T> clock_now()
	{
		return std::chrono::time_point_cast< clock_resolution_t >(T::now());
	}
	
	
	// measure duration time
	template< class T >
	constexpr clock_resolution_t clock_diff(const clock_point_t<T> & past)
	{
		return std::chrono::time_point_cast< clock_resolution_t >(clock_now<T>()) - past;
	}
	
	
	// return duration time
	template< class T, typename D >
	constexpr clock_measure_t clock_cast(const clock_point_t<T> & past)
	{
		return std::chrono::duration_cast<D>(clock_diff<T>(past)).count();
	}
	
	
	// time since epoch
	template< class T, typename D >
	constexpr clock_measure_t clock_epoch()
	{
		return std::chrono::duration_cast<D>(clock_now<T>()).clock_since_epoch();
	}
	
	
	// return as formated string
	template< class T >
	inline std::string clock_string(const clock_point_t<T> & past)
	{
		const clock_unit_u::rep u = clock_cast<T, clock_unit_u >(past);
		const clock_unit_s::rep s = clock_cast<T, clock_unit_s >(past);
		const clock_unit_m::rep m = clock_cast<T, clock_unit_m >(past);
		const clock_unit_h::rep h = clock_cast<T, clock_unit_h >(past);
		
		const std::string h_str = (h > 0) ? std::to_string(h) + "h ": std::string();
		const std::string m_str = (m > 0) ? std::to_string(m % 60) + "m ": std::string();
		const std::string s_str = (s > 0) ? std::to_string(s % 60) + "s ": std::string();
		const std::string u_str = (u > 0) ? std::to_string(u % 1000) + "ms": "<1ms";
		
		return h_str + m_str + s_str + u_str;
	}
	
	
	// return rate as formated string
	template< class T >
	inline std::string clock_rate(const clock_point_t<T> & past, const size_t & count)
	{
		return clock_string<T>(clock_diff<T>(past) / static_cast<clock_measure_t>(count));
	}
	
	
	// return expected time as formated string
	template< class T >
	inline std::string clock_expect(const clock_point_t<T> & past, const size_t & count, const size_t & full_count)
	{
		return static_cast<clock_measure_t>(full_count) / clock_rate<T>(past, count);
	}
}
*/

//
// Clock the time between initialisation and time return
//
class Clock
{
public:
	
	typedef double float_t;
	
	// add clock
	void operator << (const Clock &);
	
	// print to stream
	void print(std::ostream & = std::cout) const; // print formated string for both clocks
	
	// reset clock
	void reset(); // reset both clocks
	
	// Timestamp of current date/time
	static std::string stamp(bool = true);
	
	
private:
	
	typedef std::chrono::milliseconds tick_t; // set to milliseconds resolution
	
	typedef std::chrono::milliseconds clock_ms;
	typedef std::chrono::seconds      clock_s;
	typedef std::chrono::minutes      clock_m;
	typedef std::chrono::hours        clock_h;
	
	typedef std::chrono::duration< float_t, clock_ms::period > period_ms;
	typedef std::chrono::duration< float_t, clock_s::period  > period_s;
	typedef std::chrono::duration< float_t, clock_m::period  > period_m;
	typedef std::chrono::duration< float_t, clock_h::period  > period_h;
	
	template< class T_CLOCK >
	struct Duration
	{
		typedef std::chrono::time_point< T_CLOCK, tick_t > timepoint_t;
		
		
		// construct
		Duration()
		: timepoint(std::chrono::time_point_cast< tick_t >(T_CLOCK::now()))
		{}
		
		// return as formated string
		std::string str() const
		{
			std::ostringstream oss;
			
			const tick_t d = this->duration();
			
			const clock_ms::rep ms = std::chrono::duration_cast< clock_ms >(d).count();
			const clock_s::rep  s  = std::chrono::duration_cast< clock_s  >(d).count();
			const clock_m::rep  m  = std::chrono::duration_cast< clock_m  >(d).count();
			const clock_h::rep  h  = std::chrono::duration_cast< clock_h  >(d).count();
			
			if (h  > 0) oss << std::to_string(h)         << "h ";
			if (m  > 0) oss << std::to_string(m % 60)    << "m ";
			if (s  > 0) oss << std::to_string(s % 60)    << "s ";
			if (ms > 0) oss << std::to_string(ms % 1000) << "ms";
			else        oss << "<1ms";
			
			return(oss.str());
		}
		
		// print formated string
		void print(std::ostream & stream = std::cout) const
		{
			stream << this->get() << std::endl;
		}
		
		// reset duration
		void reset()
		{
			this->timepoint = std::chrono::time_point_cast< tick_t >(T_CLOCK::now());
		}
		
		// return durations
		//
		float_t ms() const // milliseconds
		{
			return(std::chrono::duration_cast< period_ms >(this->duration()).count());
		}
		float_t s() const  // seconds
		{
			return(std::chrono::duration_cast< period_s >(this->duration()).count());
		}
		float_t m() const  // minutes
		{
			return(std::chrono::duration_cast< period_m >(this->duration()).count());
		}
		float_t h() const  // hours
		{
			return(std::chrono::duration_cast< period_h >(this->duration()).count());
		}
		
		// measure duration time
		const tick_t duration() const
		{
			return(std::chrono::time_point_cast< tick_t >(T_CLOCK::now()) - this->timepoint);
		}
		
		// add duration
		void operator << (const Duration<T_CLOCK> & other)
		{
			this->timepoint -= other.duration();
		}
		
		timepoint_t timepoint; // start time
	};
	
	
public:
	
	Duration< std::chrono::system_clock > process; // system time
	Duration< std::chrono::steady_clock > elapsed; // user experienced time
};

#endif /* Clock_hpp */
