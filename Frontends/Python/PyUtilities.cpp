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
 * SmolDock is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SmolDock.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include "PyUtilities.h"

#include <Utilities/ReScorer.h>


namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;

/*        ReScorer(Protein &prot, Molecule &mol, std::function<double(const iConformer &, const iTransform &, const iProtein &)>& scorFunc);

        bool prepare();
        double getScore();*/

void export_Utilities() {
    p::class_<sd::ReScorer>("ReScorer", p::init<sd::Protein&, sd::Molecule&,std::function<double(const sd::iConformer &, const sd::iTransform &, const sd::iProtein &)>& >())
        .def("prepare", &sd::ReScorer::prepare)
        .def("getScore", &sd::ReScorer::getScore)
        //.def("set", &World::set)
        //.add_property("rovalue", &Num::get)
        //.add_property("value", &Num::get, &Num::set);
            ;

}

