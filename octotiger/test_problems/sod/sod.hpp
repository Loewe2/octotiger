/*
 * rotating_star.hpp
 *
 *  Created on: Oct 12, 2018
 *      Author: dmarce1
 */

#ifndef ___SOD_HPP_
#define ___SOD_HPP_

#include "octotiger/real.hpp"

#include <octotiger/debug_vector.hpp>

oct::vector<real> sod_shock_tube_init(real x0, real y, real z, real);
oct::vector<real> sod_shock_tube_analytic(real x0, real y, real z, real);


#endif /* ROTATING_STAR_ROTATING_STAR_HPP_ */
