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
#include <Engines/ConformerDockingEngine.h>

namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;

class dummyClassA {
};

p::tuple wrap_centerarray_getter(sd::Engine::AbstractDockingEngine::DockingBoxSetting &thisobj) {
    p::list a;
    a.append(thisobj.center.at(0));
    a.append(thisobj.center.at(1));
    a.append(thisobj.center.at(2));
    return p::tuple(a);
}

void wrap_centerarray_setter(sd::Engine::AbstractDockingEngine::DockingBoxSetting &thisobj, p::tuple center) {
    unsigned int tuple_len = p::extract<unsigned int>(center.attr("__len__")());
    if (tuple_len == 3) {
        thisobj.center[0] = p::extract<double>(center[0]);
        thisobj.center[1] = p::extract<double>(center[1]);
        thisobj.center[2] = p::extract<double>(center[2]);
    }
}



void export_Engines() {

    p::class_<sd::Engine::AbstractDockingEngine::DockingBoxSetting>("DockingBoxSetting")
            .def_readwrite("type", &sd::Engine::AbstractDockingEngine::DockingBoxSetting::type)
            .def_readwrite("radius", &sd::Engine::AbstractDockingEngine::DockingBoxSetting::radius)
            .add_property("center", &wrap_centerarray_getter, &wrap_centerarray_setter);

    p::enum_<sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type>("DockingBoxType")
            .value("everything", sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::everything)
                    //.value("solventExposed", sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::solventExposed)
            .value("centeredAround", sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround);

    p::enum_<sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape>("DockingBoxShape")
            .value("sphere", sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::sphere)
            .value("cube", sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::cube);

    p::enum_<sd::Score::ScoringFunctionType>("ScoringFunctionType")
            .value("VinaRigid", sd::Score::ScoringFunctionType::VinaRigid)
            .value("Vina", sd::Score::ScoringFunctionType::Vina);

    p::enum_<sd::Heuristics::GlobalHeuristicType>("GlobalHeuristicType")
            .value("RandomRestart", sd::Heuristics::GlobalHeuristicType::RandomRestart)
            .value("IteratedLocalSearch", sd::Heuristics::GlobalHeuristicType::IteratedLocalSearch)
            .value("OnlyLocal", sd::Heuristics::GlobalHeuristicType::OnlyLocal)
            .value("SimulatedAnnealing", sd::Heuristics::GlobalHeuristicType::SimulatedAnnealing);

    p::enum_<sd::Optimizer::LocalOptimizerType>("LocalOptimizerType")
            .value("L_BFGS", sd::Optimizer::LocalOptimizerType::L_BFGS)
            .value("GradientDescentLineSearch", sd::Optimizer::LocalOptimizerType::GradientDescentLineSearch);

    {
        p::scope engineScope
                = p::class_<dummyClassA>("Engine");

        /*             // /// Actions /////////////
        bool setupDockingEngine() final;

        void runDockingEngine() final;*/
        p::class_<sd::Engine::ConformerDockingEngine>("ConformerDockingEngine", p::init<
                unsigned int,
                unsigned int,
                sd::Protein *,
                sd::Molecule *,
                sd::Score::ScoringFunctionType,
                sd::Heuristics::GlobalHeuristicType,
                sd::Optimizer::LocalOptimizerType,
                unsigned int
        >())
                .def("setupDockingEngine", &sd::Engine::ConformerDockingEngine::setupDockingEngine)
                .def("runDockingEngine", &sd::Engine::ConformerDockingEngine::runDockingEngine)
                .def("getDockingResult", &sd::Engine::ConformerDockingEngine::getDockingResult)
                .def("setDockingBox", &sd::Engine::ConformerDockingEngine::setDockingBox)
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