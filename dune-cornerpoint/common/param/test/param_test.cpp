//===========================================================================
//
// File: param_test.cpp
//
// Created: Sun Dec 13 20:08:36 2009
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


#define BOOST_TEST_DYN_LINK
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE ParameterTest
#include <boost/test/unit_test.hpp>

#include "../ParameterGroup.hpp"

using namespace Dune;

BOOST_AUTO_TEST_CASE(commandline_syntax_init)
{
    typedef char* cp;
    const int argc = 8;
    cp argv[argc] = { "program_command",
                      "topitem=somestring",
                      "/slashtopitem=anotherstring",
                      "/group/item=1",
                      "/group/anotheritem=2",
                      "/group/subgroup/item=3",
                      "/group/subgroup/anotheritem=4",
                      "/group/item=overridingstring" };
    parameter::ParameterGroup p(argc, argv);
    BOOST_CHECK(p.get<std::string>("topitem") == "somestring");

    // Tests that only run in debug mode.
#ifndef NDEBUG
#endif
}
