//===========================================================================
//
// File: Average.hpp<2>
//
// Created: Wed Jul  1 12:25:57 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009 SINTEF ICT, Applied Mathematics.
  Copyright 2009 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_AVERAGE_HEADER
#define OPENRS_AVERAGE_HEADER


#include <boost/static_assert.hpp>
#include <tr1/type_traits>
#include <cmath>

namespace Dune {
    namespace utils {



	/// Computes the arithmetic average:
	/// a_A(t_1, t_2) = (t_1 + t_2)/2.
	template <typename T, typename Tresult = T>
	Tresult arithmeticAverage(const T& t1, const T& t2)
	{
	    // To avoid some user errors, we disallow taking averages of
	    // integral values.
	    BOOST_STATIC_ASSERT(std::tr1::is_integral<T>::value == false);
	    Tresult retval(t1);
	    retval += t2;
	    retval *= 0.5;
	    return retval;
	}



	/// Computes the geometric average:
	/// a_G(t_1, t_2) = \sqrt{t_1 t_2}.
	/// Since we use std::sqrt(), this function will only
	/// compile when T is convertible to and from double.
	template <typename T>
	T geometricAverage(const T& t1, const T& t2)
	{
	    // To avoid some user errors, we disallow taking averages of
	    // integral values.
	    BOOST_STATIC_ASSERT(std::tr1::is_integral<T>::value == false);
	    return std::sqrt(t1*t2);
	}



	/// Computes the harmonic average:
	/// a_H(t_1, t_2) = \frac{2}{1/t_1 + 1/t_2} = \frac{2 t_1 t_2}{t_1 + t_2}.
	template <typename T>
	T harmonicAverage(const T& t1, const T& t2)
	{
	    // To avoid some user errors, we disallow taking averages of
	    // integral values.
	    BOOST_STATIC_ASSERT(std::tr1::is_integral<T>::value == false);
	    return (2*t1*t2)/(t1 + t2);
	}



    } // namespace utils
} // namespace Dune



#endif // OPENRS_AVERAGE_HEADER
