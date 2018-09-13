//
//  Copyright (c) 2018 Patrick K. Albers. All rights reserved.
//

#ifndef Decimal_h
#define Decimal_h

#include <limits>
#include <vector>


using decimal_t = double;
using decimal_vector_t = std::vector< decimal_t >;

static constexpr decimal_t decimal_eps = std::numeric_limits<decimal_t>::epsilon();
static constexpr decimal_t decimal_min = std::numeric_limits<decimal_t>::min();
static constexpr decimal_t decimal_max = std::numeric_limits<decimal_t>::max();

static constexpr decimal_t decimal_neg = std::numeric_limits<decimal_t>::lowest();
static constexpr decimal_t decimal_pos = std::numeric_limits<decimal_t>::max();

static constexpr decimal_t decimal_nil = static_cast<decimal_t>(0);
static constexpr decimal_t decimal_one = static_cast<decimal_t>(1);
static constexpr decimal_t decimal_two = static_cast<decimal_t>(2);
static constexpr decimal_t decimal_six = static_cast<decimal_t>(6);
static constexpr decimal_t decimal_100 = static_cast<decimal_t>(100);

static constexpr decimal_t decimal_err = static_cast<decimal_t>(1e-08);


#endif /* Decimal_h */

