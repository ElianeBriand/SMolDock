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

#include <iostream>
#include <memory>


#include "Structures/Molecule.h"
#include "Structures/Protein.h"
#include "Structures/Atom.h"
#include "Utilities/DockingResultPrinter.h"
#include "Engines/ConformerDockingEngine.h"
#include "Engines/ScoringFunctions/VinaLikeRigid.h"
#include <Engines/Internals/InternalsUtilityFunctions.h>

#include <Utilities/PDBWriter.h>
#include <Utilities/CSVReader.h>
#include <Utilities/SMARTSMatcher.h>

#include <Structures/InputModifiers/InputModifierInterface.h>
#include <Structures/InputModifiers/VinaCompatibility.h>

#include <Utilities/Version.h>
#include <Utilities/Calibration/Calibrator.h>

#include <thread>

#include <boost/lexical_cast.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>

#define USE_BOOST_KARMA

#include <bprinter/table_printer.h>

#include <tbb/task_scheduler_init.h>

namespace sd = SmolDock;

#include <boost/mpi.hpp>
#include <iostream>
#include <string>
#include <boost/serialization/string.hpp>

namespace mpi = boost::mpi;

#include <Frontends/FrontendCommon.h>
#include <Utilities/Calibration/MPICalibratorDirector.h>
#include <Utilities/Calibration/MPICalibratorNode.h>

#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

int main() {


    //        tbb::task_scheduler_init tbbInit(std::thread::hardware_concurrency());


    mpi::environment env;
    mpi::communicator world;

    if (world.rank() == 0) {
        sd::setupLogPrinting(false, false, "[MPI_0] ");
        auto cdirector = std::make_shared<sd::Calibration::MPICalibratorDirector>(env, world,
                                                                                  sd::Score::ScoringFunctionType::VinaCovalentReversible,
                                                                                  sd::Heuristics::GlobalHeuristicType::SimulatedAnnealing,
                                                                                  sd::Optimizer::LocalOptimizerType::L_BFGS,
                                                                                  1000, // Max epoch
                                                                                  0.5, // Initial learning rate
                                                                                  32574, // RNG seed
                                                                                  5, // num generated starting conformer per ligand
                                                                                  3, // num retry per conformer (best score is kept)
                                                                                  18, //batch size
                                                                                  sd::Heuristics::emptyParameters);

        //cdirector->coefficientsToCalibrate({"Gauss1","Gauss2","RepulsionExceptCovalent","Hydrophobic","Hydrogen","CovalentReversible"});
        cdirector->coefficientsToCalibrate({"CovalentReversible"});

        sd::Engine::AbstractDockingEngine::DockingBoxSetting setting;
        setting.type = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround;
        setting.center = {2.1, -5.6, 52.0};
        setting.radius = 14.0;
        sd::Calibration::Calibrator::ReceptorID recID1 =
                cdirector->addReceptorFromFile("/home/briand/CLionProjects/SmolDock/DockingTests/hCES1/1YA8.pdb", setting);

        cdirector->applySpecialResidueTypingFromRecID(recID1,
                sd::AminoAcid::AAType::serine,
                221,
                sd::SpecialResidueTyping::covalentReversibleSerineOH);

        sd::CSVReader chembl_csv("/home/briand/CLionProjects/SmolDock/DockingTests/hCES1/chembl_data.tsv","\t",true);
        std::vector<std::map<std::string,std::string>> chembl_data = chembl_csv.getRowsAsMap();
        sd::SMARTSMatcher matchDiKetone("[c,C]C(=O)C(=O)[c,C]");

        for (unsigned int j = 0; j < chembl_data.size(); ++j) {
            std::string smiles = chembl_data.at(j).at("CANONICAL_SMILES");
            double ki = boost::lexical_cast<double>(chembl_data.at(j).at("STANDARD_VALUE")) * 1e-9;

            if(matchDiKetone.matchesSMILES(smiles) && ki < 100e-9)
                cdirector->addReferenceLigand_SMILES_Ki(recID1,smiles, ki);
        }

        cdirector->applyVariantToAllLigands("[c,C][C:1](=O)[c,C]", sd::Atom::AtomVariant::covalentReversibleAcceptor);

        cdirector->setupCalibration();
        cdirector->runCalibration();

    } else {
        tbb::task_scheduler_init tbbInit(std::thread::hardware_concurrency());
        //tbb::task_scheduler_init tbbInit(1);
        sd::setupLogPrinting(false, false, "[MPI_" + lexical_cast<std::string>(world.rank()) +"] ");
        auto cnode = std::make_shared<sd::Calibration::MPICalibratorNode>(env, world);
        cnode->runNode();
    }


    return 0;


}
