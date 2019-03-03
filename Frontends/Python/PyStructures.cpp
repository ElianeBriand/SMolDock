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

#include <memory>

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Structures/Molecule.h>
#include <Structures/Protein.h>
#include <Structures/Results/DockingResult.h>

namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;

// Class: Molecule
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(populateFromPDBFile_overloads, sd::Molecule::populateFromPDBFile, 1, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(populateFromSMILES_overloads, sd::Molecule::populateFromSMILES, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(generateConformer_overloads, sd::Molecule::generateConformer, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(populateFromMolFile_overloads, sd::Molecule::populateFromMolFile, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(populateFromMol2File_overloads, sd::Molecule::populateFromMol2File, 1, 3)

// Class: Protein
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(prot_populateFromPDB_overloads, sd::Protein::populateFromPDB, 1, 2)


// Remember : Any function added to a class whose initial argument matches the class (or any base) will act like a member function in Python.
// Aka defining f(MyClass*, arg1, arg2,) allow you to easily patch-in python-specific member function without altering the C++ class (just adding friend functions)

namespace SmolDock::Wrapper {

    RDKit::ROMol *pywrap_returnROMolFromMol(const SmolDock::Molecule &mol) {
        auto ret = new RDKit::ROMol(static_cast<RDKit::ROMol>(*mol.rwmol));
        return ret;
    }

}

/*        bool populateFromMolFile(const std::string &filename, unsigned int seed = 36754,
                                 std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});

        bool populateFromMol2File(const std::string &filename, unsigned int seed = 36754,
                                  std::vector<std::shared_ptr<InputPostProcessor::InputPostProcessor> > postProcessors = {});
*/

void export_Structures() {
    p::class_<sd::Molecule>("Molecule")
            .def("populateFromPDB", &sd::Molecule::populateFromPDBFile, populateFromPDBFile_overloads())
            .def("populateFromSMILES", &sd::Molecule::populateFromSMILES, populateFromSMILES_overloads())
            .def("populateFromMol", &sd::Molecule::populateFromMolFile, populateFromMolFile_overloads())
            .def("populateFromMol2", &sd::Molecule::populateFromMol2File, populateFromMol2File_overloads())
            .add_property("numberOfAtoms", &sd::Molecule::numberOfAtoms)
            .add_property("numberOfBonds", &sd::Molecule::numberOfBonds)
            .def("getInitialConformer", &sd::Molecule::getInitialConformer)
            .def("generateConformer", &sd::Molecule::generateConformer, generateConformer_overloads())
            .add_property("residueName", &sd::Molecule::getResidueName, &sd::Molecule::setResidueName)
            .def("updateAtomPositionsFromiConformer", &sd::Molecule::updateAtomPositionsFromiConformer)
            .def("deepcopy", &sd::Molecule::deepcopy)
            .def("getRDKitMol", &sd::Wrapper::pywrap_returnROMolFromMol, p::return_value_policy<p::manage_new_object>())
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
            .def_readwrite("variant", &sd::iConformer::variant);

    p::class_<sd::iProtein>("iProtein")
            .def_readwrite("center_x", &sd::iProtein::center_x)
            .def_readwrite("center_y", &sd::iProtein::center_y)
            .def_readwrite("center_z", &sd::iProtein::center_z)
            .def_readwrite("x", &sd::iProtein::x)
            .def_readwrite("y", &sd::iProtein::y)
            .def_readwrite("z", &sd::iProtein::z)
            .def_readwrite("type", &sd::iProtein::type)
            .def_readwrite("variant", &sd::iProtein::variant)
            .def_readwrite("AAId_to_AtomPositionInVect", &sd::iProtein::AAId_to_AtomPositionInVect);


    p::class_<sd::Protein>("Protein", p::init<>())
            .def("populateFromPDB", &sd::Protein::populateFromPDB, prot_populateFromPDB_overloads())
            .def("getiProtein", &sd::Protein::getiProtein);

    p::class_<sd::DockingResult, std::shared_ptr<sd::DockingResult> >("DockingResult")
            .def_readwrite("ligandPoses", &sd::DockingResult::ligandPoses);


}
