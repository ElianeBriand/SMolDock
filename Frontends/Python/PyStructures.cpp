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

#include <Structures/Molecule.h>
#include <Structures/Protein.h>
#include <Structures/Results/DockingResult.h>

namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(populateFromPDB_overloads, sd::Molecule::populateFromPDB, 1, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(populateFromSMILES_overloads, sd::Molecule::populateFromSMILES, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generateConformer_overloads, sd::Molecule::generateConformer, 1, 3)



void export_Structures()
{
    p::class_<sd::Molecule>("Molecule")
            .def("populateFromPDB", &sd::Molecule::populateFromPDB, populateFromPDB_overloads())
            .def("populateFromSMILES", &sd::Molecule::populateFromSMILES, populateFromSMILES_overloads())
            .add_property("numberOfAtoms", &sd::Molecule::numberOfAtoms)
            .add_property("numberOfBonds", &sd::Molecule::numberOfBonds)
            .def("getInitialConformer", &sd::Molecule::getInitialConformer)
            .def("generateConformer", &sd::Molecule::generateConformer, generateConformer_overloads() )
            .add_property("residueName", &sd::Molecule::getResidueName, &sd::Molecule::setResidueName)
            .def("updateAtomPositionsFromiConformer", &sd::Molecule::updateAtomPositionsFromiConformer)
            .def("deepcopy", &sd::Molecule::deepcopy)
            //.def("populateFromPDB", &sd::Molecule::populateFromPDB)
            //.def("set", &World::set)
            //.add_property("rovalue", &Num::get)
            //.add_property("value", &Num::get, &Num::set);
            ;

    p::class_<sd::iConformer>("iConformer")
            .def_readwrite("x", &sd::iConformer::x)
            .def_readwrite("y", &sd::iConformer::y)
            .def_readwrite("z", &sd::iConformer::z)
            .def_readwrite("type", &sd::iConformer::type)
            .def_readwrite("variant", &sd::iConformer::variant)
            ;

    p::class_<sd::iProtein>("iProtein")
            .def_readwrite("center_x", &sd::iProtein::center_x)
            .def_readwrite("center_y", &sd::iProtein::center_y)
            .def_readwrite("center_z", &sd::iProtein::center_z)
            .def_readwrite("x", &sd::iProtein::x)
            .def_readwrite("y", &sd::iProtein::y)
            .def_readwrite("z", &sd::iProtein::z)
            .def_readwrite("type", &sd::iProtein::type)
            .def_readwrite("variant", &sd::iProtein::variant)
            .def_readwrite("AAId_to_AtomPositionInVect", &sd::iProtein::AAId_to_AtomPositionInVect)
            ;


    p::class_<sd::Protein>("Protein")
            .def("populateFromPDB", &sd::Protein::populateFromPDB)
            .def("getiProtein", &sd::Protein::getiProtein)
            ;

    p::class_<sd::DockingResult>("DockingResult")
            .def_readwrite("ligandPoses", &sd::DockingResult::ligandPoses)
            ;


}
