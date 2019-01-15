//
// Created by eliane on 11/01/19.
//

#include "PyEngine.h"


#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include <Structures/Molecule.h>
#include <Structures/Protein.h>
#include <Structures/Results/DockingResult.h>

#include <Engines/ScoringFunctions/ScoringFunctionFactory.h>
#include <Engines/GlobalHeuristics/HeuristicFactory.h>
#include <Engines/LocalOptimizers/OptimizerFactory.h>
#include <Engines/ConformerRigidDockingEngine.h>

namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;

class dummyClassA {
};


void export_Engines() {
    p::enum_<sd::Score::ScoringFunctionType>("ScoringFunctionType")
            .value("VinaRigid", sd::Score::ScoringFunctionType::VinaRigid);

    p::enum_<sd::Heuristics::GlobalHeuristicType>("GlobalHeuristicType")
            .value("RandomRestart", sd::Heuristics::GlobalHeuristicType::RandomRestart);

    p::enum_<sd::Optimizer::LocalOptimizerType>("LocalOptimizerType")
            .value("L_BFGS", sd::Optimizer::LocalOptimizerType::L_BFGS)
            .value("GradientDescentLineSearch", sd::Optimizer::LocalOptimizerType::GradientDescentLineSearch);

    {
        p::scope engineScope
                = p::class_<dummyClassA>("Engine");

        /*             // /// Actions /////////////
        bool setupDockingEngine() final;

        void runDockingEngine() final;*/
        p::class_<sd::Engine::ConformerRigidDockingEngine>("ConformerRigidDockingEngine", p::init<
                unsigned int,
                sd::Protein*,
                sd::Molecule*,
                sd::Score::ScoringFunctionType,
                sd::Heuristics::GlobalHeuristicType,
                sd::Optimizer::LocalOptimizerType,
                unsigned int
        >())
                .def("setupDockingEngine", &sd::Engine::ConformerRigidDockingEngine::setupDockingEngine)
                .def("runDockingEngine", &sd::Engine::ConformerRigidDockingEngine::runDockingEngine)
                .def("getDockingResult", &sd::Engine::ConformerRigidDockingEngine::getDockingResult)
            /*   .def("populateFromPDB", &sd::Molecule::populateFromPDB, populateFromPDB_overloads())
               .def("populateFromSMILES", &sd::Molecule::populateFromSMILES, populateFromSMILES_overloads())
               .add_property("numberOfAtoms", &sd::Molecule::numberOfAtoms)
               .add_property("numberOfBonds", &sd::Molecule::numberOfBonds)
               .def("getInitialConformer", &sd::Molecule::getInitialConformer)
               .def("generateConformer", &sd::Molecule::generateConformer, generateConformer_overloads() )
               .add_property("residueName", &sd::Molecule::getResidueName, &sd::Molecule::setResidueName)
               .def("updateAtomPositionsFromiConformer", &sd::Molecule::updateAtomPositionsFromiConformer)
               .def("deepcopy", &sd::Molecule::deepcopy)*/
            //.def("populateFromPDB", &sd::Molecule::populateFromPDB)
            //.def("set", &World::set)
            //.add_property("rovalue", &Num::get)
            //.add_property("value", &Num::get, &Num::set);
                ;


    }


}