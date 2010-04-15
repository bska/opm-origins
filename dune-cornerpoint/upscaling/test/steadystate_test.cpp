//===========================================================================
//
// File: steadystate_test.cpp
//
// Created: Fri Aug 28 14:11:03 2009
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
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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


//#define VERBOSE

#include <dune/upscaling/SteadyStateUpscaler.hpp>
#include <dune/common/Units.hpp>
#include <dune/common/SparseTable.hpp>
#include <cmath>

using namespace Dune;
using namespace Dune::prefix;
using namespace Dune::unit;


/// Reads saturation and pressure drop data from an input stream.
template <class Istream>
void readControl(Istream& is, std::vector<double>& saturations, SparseTable<double>& all_pdrops)
{
    int num_sats;
    is >> num_sats;
    std::vector<double> sat(num_sats);
    SparseTable<double> ap;
    std::vector<double> tmp_pd;
    for (int i = 0; i < num_sats; ++i) {
        is >> sat[i];
        int num_pdrops;
        is >> num_pdrops;
        tmp_pd.resize(num_pdrops);
        for (int j = 0; j < num_pdrops; ++j) {
            is >> tmp_pd[j];
        }
        ap.appendRow(tmp_pd.begin(), tmp_pd.end());
    }
    all_pdrops.swap(ap);
    saturations.swap(sat);
}



/// Writes saturation and pressure drop data to an output stream.
template <class Ostream>
void writeControl(Ostream& os, const std::vector<double>& saturations, const SparseTable<double>& all_pdrops)
{
    int num_sats = saturations.size();
    os << num_sats << '\n';
    for (int i = 0; i < num_sats; ++i) {
        os << saturations[i] << '\n';
        int num_pdrops = all_pdrops[i].size();
        os << num_pdrops;
        for (int j = 0; j < num_pdrops; ++j) {
            os << ' ' << all_pdrops[i][j];
        }
        os << '\n';
    }
}



int main(int argc, char** argv)
{
    // Initialize.
    parameter::ParameterGroup param(argc, argv);
    // MPIHelper::instance(argc,argv);

    // Control structure.
    std::vector<double> saturations;
    SparseTable<double> all_pdrops;
    bool from_file = param.has("sat_pdrop_filename");
    if (from_file) {
        std::string filename = param.get<std::string>("sat_pdrop_filename");
        std::ifstream file(filename.c_str());
        if (!file) {
            THROW("Could not open file " << filename);
        }
        readControl(file, saturations, all_pdrops);
    } else {
        // Get a linear range of saturations.
        int num_sats = param.getDefault("num_sats", 4);
        double min_sat = param.getDefault("min_sat", 0.2);
        double max_sat = param.getDefault("max_sat", 0.8);
        saturations.resize(num_sats);
        for (int i = 0; i < num_sats; ++i) {
            double factor = num_sats == 1 ? 0 : double(i)/double(num_sats - 1);
            saturations[i] = (1.0 - factor)*min_sat + factor*max_sat;
        }
        // Get a logarithmic range of pressure drops.
        int num_pdrops = param.getDefault("num_pdrops", 5);
        double log_min_pdrop = std::log(param.getDefault("min_pdrop", 1e2));
        double log_max_pdrop = std::log(param.getDefault("max_pdrop", 1e6));
        std::vector<double> pdrops;
        pdrops.resize(num_pdrops);
        for (int i = 0; i < num_pdrops; ++i) {
            double factor = num_pdrops == 1 ? 0 : double(i)/double(num_pdrops - 1);
            pdrops[i] = std::exp((1.0 - factor)*log_min_pdrop + factor*log_max_pdrop);
        }
        // Assign the same pressure drops to all saturations.
        for (int i = 0; i < num_sats; ++i) {
            all_pdrops.appendRow(pdrops.begin(), pdrops.end());
        }
    }

    // Print the saturations and pressure drops.
    // writeControl(std::cout, saturations, all_pdrops);

    // Initialize upscaler.
    SteadyStateUpscaler upscaler;
    upscaler.init(param);

    // First, compute an upscaled permeability.
    SteadyStateUpscaler::permtensor_t upscaled_K = upscaler.upscaleSinglePhase();
    SteadyStateUpscaler::permtensor_t upscaled_K_copy = upscaled_K;
    upscaled_K_copy *= (1.0/(milli*darcy));
    std::cout.precision(15);
    std::cout << "Upscaled K in millidarcy:\n" << upscaled_K_copy << std::endl;

    // Then, compute some upscaled relative permeabilities.
    int num_cells = upscaler.grid().size(0);
    /*
//     const int num_sats = 1;
//     double saturations[num_sats] = { 0.3 };
    const int num_sats = 7;
    double saturations[num_sats] = { 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
    // const int num_sats = 11;
    // double saturations[num_sats] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    const int num_pdrops = 5;
    double pdrops[num_pdrops] = { 1e1, 1e2, 1e3, 1e4, 1e5};
//     const int num_pdrops = 1;
//     double pdrops[num_pdrops] = { 1e1 };
*/
    int num_sats = saturations.size();
    for (int i = 0; i < num_sats; ++i) {
	// Starting every computation with a trio of uniform profiles.
	std::vector<double> init_sat(num_cells, saturations[i]);
	boost::array<std::vector<double>, 3> init_sats = {{ init_sat, init_sat, init_sat }};
        const SparseTable<double>::row_type pdrops = all_pdrops[i];
        int num_pdrops = pdrops.size();
	for (int j = 0; j < num_pdrops; ++j) {
	    double pdrop = pdrops[j];
	    SteadyStateUpscaler::permtensor_t upscaled_relperm
		= upscaler.upscaleSteadyState(init_sats, saturations[i], pdrop, upscaled_K);
            boost::array<double, 3> usats = upscaler.lastSaturationsUpscaled();
	    std::cout << "\n\nTensor of upscaled relperms for initial saturation " << saturations[i]
                      << ", real steady-state saturations " << usats[0] << ' ' << usats[1] << ' ' << usats[2]
                      << " and pressure drop " << pdrop << ":\n" << upscaled_relperm << "\n" << std::endl;
	    // Changing initial saturations for next pressure drop to equal the steady state of the last
	    init_sats = upscaler.lastSaturations();
	}
    }
}
