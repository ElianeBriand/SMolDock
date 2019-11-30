//
// Created by eliane on 15/03/19.
//

#include "PyCalibrator.h"


#include <boost/python.hpp>

#include <Structures/Molecule.h>
#include <Structures/Protein.h>
#include <Structures/Results/DockingResult.h>

#include <Engines/ScoringFunctions/ScoringFunctionFactory.h>
#include <Engines/GlobalHeuristics/HeuristicFactory.h>
#include <Engines/LocalOptimizers/OptimizerFactory.h>
#include <Engines/ConformerDockingEngine.h>

#include <Utilities/Calibration/Calibrator.h>

namespace p = boost::python;

namespace sd = SmolDock;


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(addReferenceLigand_SMILES_Ki_overloads, sd::Calibration::Calibrator::addReferenceLigand_SMILES_Ki, 3, 4)


void export_Calibrator() {

    //p::class_<sd::Calibration::Calibrator::ReceptorID>("ReceptorID");

    p::class_<sd::Calibration::Calibrator>("Calibrator", p::init<
            sd::Score::ScoringFunctionType,
            sd::Heuristics::GlobalHeuristicType,
            sd::Optimizer::LocalOptimizerType,
            p::optional<unsigned int,
                    double,
                    unsigned int,
                    unsigned int,
                    unsigned int,
                    unsigned int,
                    sd::Heuristics::HeuristicParameters> >())
            .def("addReceptor", &sd::Calibration::Calibrator::addReceptor)
            .def("addReferenceLigand_SMILES_Ki", &sd::Calibration::Calibrator::addReferenceLigand_SMILES_Ki, addReferenceLigand_SMILES_Ki_overloads())
//            .def("runDockingEngine", &sd::Engine::ConformerDockingEngine::runDockingEngine)
//            .def("getDockingResult", &sd::Engine::ConformerDockingEngine::getDockingResult)
//            .def("setDockingBox", &sd::Engine::ConformerDockingEngine::setDockingBox)
            ;


}