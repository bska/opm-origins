//===========================================================================
//
// File: SpecialEclipseFields.hpp
//
// Created: Mon Sep 21 14:09:54 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_SPECIALECLIPSEFIELDS_HEADER
#define OPENRS_SPECIALECLIPSEFIELDS_HEADER

#include <string>
#include <fstream>
#include <limits>
#include <dune/common/ErrorMacros.hpp>
#include <iterator>
#include "EclipseGridParserHelpers.hpp"

namespace Dune
{


// Abstract base class for special fields. A special field is a field
// with more than one data type.
struct SpecialBase {
    virtual ~SpecialBase() {}                       // Default destructor
    //virtual std::string name() const = 0;           // Keyword name
    virtual void read(std::istream& is) = 0;        // Reads data
    //virtual void write(std::ostream& os) const = 0; // Writes data
};




/// Class for keyword SPECGRID 
struct SPECGRID : public SpecialBase
{
    std::vector<int> dimensions; // Number of grid blocks in x-, y- and z-directions.
    int numres;          // Number of reservoirs. 
    char qrdial;         // Coordinates. F=cartesian, T=Cylindrical(radial).

    SPECGRID()
    {
	dimensions.resize(3,1);
	numres = 1;
	qrdial = 'F';
    }

    virtual ~SPECGRID()
    {
    }

    virtual std::string name() const {return std::string("SPECGRID");}

    virtual void read(std::istream& is)
    {
	const int ndim = 3;
	std::vector<int> data(ndim+1,1);
	int nread = readDefaultedVectorData(is , data, ndim+1);
	int nd = std::min(nread, ndim);
	copy(data.begin(), data.begin()+nd, &dimensions[0]);
	numres = data[ndim];
	std::string candidate;
	is >> candidate;
	if (candidate == "/") {
	    return;
	} else {
	    qrdial = candidate[0];
	}
	
	if (ignoreSlashLine(is)) {
	    return;
	} else {
	    THROW("End of file reading" << name());
	}
    }
	
    virtual void write(std::ostream& os) const
    {
	os << name() << std::endl;
	os << dimensions[0] << " " << dimensions[1]  << " "
	   << dimensions[2] << " " << numres << " " << qrdial << std::endl;
	os << std::endl;
    }
};




/// Class holding segment data of keyword FAULTS.
struct FaultSegment
{
    std::string fault_name;          // Fault name
    std::vector<int> ijk_coord;      // ijk-coordinates of segment cells
    std::string face;                // Fault face of cells
};

/// Class for keyword FAULTS.
struct FAULTS : public SpecialBase
{
    std::vector<FaultSegment> faults;

    FAULTS()
    {
    }

    virtual ~FAULTS()
    {
    }

    virtual std::string name() const {return std::string("FAULTS");}

    virtual void read(std::istream& is)
    {
	while(is) {
	    std::string name;
	    is >> name;
	    if (name[0] == '/') {
		break;
	    }
	    while (name.find("--") == 0) {
		// This line is a comment
		is >> ignoreLine >> name;
	    }
	    FaultSegment fault_segment;
	    fault_segment.ijk_coord.resize(6);
	    fault_segment.fault_name = name;
	    int nread = readDefaultedVectorData(is, fault_segment.ijk_coord, 6);
	    if (nread != 6) {
		THROW("Error reading fault_segment " << name);
	    }
	    is >> fault_segment.face;
	    faults.push_back(fault_segment);
	    ignoreSlashLine(is);
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << std::endl;
	for (int i=0; i<(int)faults.size(); ++i) {
	    os << faults[i].fault_name << "  ";
	    copy(faults[i].ijk_coord.begin(), faults[i].ijk_coord.end(),
		 std::ostream_iterator<int>(os, " "));
	    os << faults[i].face << std::endl;
	}
	os << std::endl;
    }
};




/// Class holding a data line of keyword MULTFLT
struct MultfltLine
{
    std::string fault_name;         // Fault name, as in FAULTS
    double transmis_multiplier;     // Transmissibility multiplier
    double diffusivity_multiplier;  // Diffusivity multiplier;
    MultfltLine() :
	fault_name(""), transmis_multiplier(1.0), diffusivity_multiplier(1.0)
    {
    }
};

/// Class for keyword MULFLT
struct MULTFLT : public SpecialBase
{
    std::vector<MultfltLine> multflts;

    MULTFLT()
    {
    }

    virtual ~MULTFLT()
    {
    }

    virtual std::string name() const {return std::string("MULTFLT");}

    virtual void read(std::istream& is)
    {
	while(is) {
	    std::string name;
	    is >> name;
	    if (name[0] == '/') {
		break;
	    }
	    while (name == "--") {
		// This line is a comment
		is >> ignoreLine >> name;
	    }
	    MultfltLine multflt_line;
	    multflt_line.fault_name = name;
	    std::vector<double> data(2,1.0);
	    if (readDefaultedVectorData(is, data, 2) == 2) {
		ignoreSlashLine(is);
	    }
	    multflt_line.transmis_multiplier = data[0];
	    multflt_line.diffusivity_multiplier = data[1];
	    multflts.push_back(multflt_line);
	}
    }

    virtual void write(std::ostream& os) const
    {
	os << name() << std::endl;
	for (int i=0; i<(int)multflts.size(); ++i) {
	    os << multflts[i].fault_name << "  " 
	       << multflts[i].transmis_multiplier << "  "
	       << multflts[i].diffusivity_multiplier <<	std::endl;
	}
	os << std::endl;
    }
};




struct TITLE : public SpecialBase
{
    std::string title;
    virtual std::string name() const
    { return std::string("TITLE"); }
    virtual void read(std::istream& is)
    { is >> ignoreLine; std::getline(is, title); }
    virtual void write(std::ostream& os) const
    { os << name() << '\n' << title << '\n'; }

};




struct START : public SpecialBase
{
    boost::gregorian::date date;
    virtual std::string name() const
    { return std::string("START"); }
    virtual void read(std::istream& is)
    { date = readDate(is); }
    virtual void write(std::ostream& os) const
    { os << name() << '\n' << date << '\n'; }
};




struct DATES : public SpecialBase
{
    std::vector<boost::gregorian::date> dates;
    virtual std::string name() const
    { return std::string("DATES"); }
    virtual void read(std::istream& is)
    {
	while(is) {
	    dates.push_back(readDate(is));
	    is >> ignoreWhitespace;
	    if (is.peek() == int('/')) {
		is >> ignoreLine;
		break;
	    }
	}
    }
    virtual void write(std::ostream& os) const
    {
	os << name() << '\n';
	copy(dates.begin(), dates.end(),
	     std::ostream_iterator<boost::gregorian::date>(os, "\n"));
    }

};




struct MultRec : public SpecialBase
{
    virtual void read(std::istream& is)
    {
#ifdef VERBOSE
        std::cout << "(dummy implementation)" << std::endl;
#endif
	const std::ctype<char>& ct = std::use_facet< std::ctype<char> >(std::locale::classic());
        while (true) {
            is >> ignoreSlashLine;
            is >> ignoreWhitespace;
            char c;
            is.get(c);
            is.putback(c);
            if (ct.is(std::ctype_base::alpha, c)) {
                break;
            }
        }
    }
};


struct SWFN : public MultRec {};
struct SOF2 : public MultRec {};
struct EQUIL : public MultRec {};
struct WELSPECS : public MultRec {};
struct COMPDAT : public MultRec {};
struct WCONINJE : public MultRec {};
struct TUNING : public MultRec {};


} // End of namespace Dune

#endif // OPENRS_SPECIALECLIPSEFIELDS_HEADER

