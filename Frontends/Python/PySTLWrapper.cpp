/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "PySTLWrapper.h"


#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <boost/python/suite/indexing/indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


#include <Structures/Molecule.h>


namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;



void export_STLWrapper()
{
    p::class_<std::vector<double> >("DoubleVect")
            .def(p::vector_indexing_suite<std::vector<double> >());

    p::class_<std::vector<sd::Molecule> >("MoleculeVect")
            .def(p::vector_indexing_suite<std::vector<sd::Molecule> >());
}