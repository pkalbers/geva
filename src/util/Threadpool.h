//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Threadpool_h
#define Threadpool_h

#include <functional>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <vector>

#include "Progress.hpp"
#include "Clock.hpp"


template<class M>
class Threadpool
{
public:
	
	using method_t = M;
	
	
	// construct
	Threadpool(const size_t & n, void (M::*f)())
	: core(n)
	, core_size(n)
	, pool_size(0)
	, func_exec(f)
	, error(nullptr)
	{}
	
	// append instance
	void task(method_t && instance)
	{
		std::lock_guard<std::mutex> lock(this->guard);
		
		this->pool.push(std::move(instance));
		++this->pool_size;
	}
	
	// run over all instances
	void exec(Progress * prog = nullptr, Clock * time = nullptr)
	{
		Clock local;
		
		try
		{
			while (!this->error)
			{
				this->guard.lock();
				
				if (this->pool_size > 0)
				{
					method_t instance = std::move(this->pool.front()); // move instance into scope
					
					this->pool.pop(); // remove from pool
					
					--this->pool_size;
					
					if (prog)
						prog->update();
					
					this->guard.unlock();
					
					(instance.*func_exec)();
					
					continue;
				}
				
				this->guard.unlock(); // unlock when empty
				break;
			}
		}
		catch (const std::exception & ex)
		{
			this->error = std::current_exception();
		}
		
		
		std::lock_guard<std::mutex> lock(this->guard);
		
		if (time)
			*time << local;
	}
	
	// run N threads
	void open(Progress * prog = nullptr, Clock * time = nullptr)
	{
		for (size_t i = 0; i < this->core_size; ++i)
		{
			this->core[i] = std::thread(&Threadpool<M>::exec, this, prog, time);
		}
	}
	
	// join threads
	void wait()
	{
		for (size_t i = 0; i < this->core_size; ++i)
		{
			this->core[i].join();
		}
		
		if (this->error)
		{
			std::rethrow_exception(this->error);
		}
	}
	
	
private:
	
	using thread_core_t = std::vector< std::thread >;
	using thread_pool_t = std::queue< method_t >;
	
	typedef void (M::*execute_f)();
	
	thread_core_t core; // threads
	thread_pool_t pool; // pool of method instances
	
	size_t core_size; // number of threads
	size_t pool_size; // number of tasks
	
	execute_f func_exec;
	
	std::exception_ptr error; // exception pointer
	std::mutex         guard; // mutex
};


#endif /* Threadpool_h */

