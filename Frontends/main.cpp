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
#include "Utilities/DockingResultPrinter.h"
#include "Engines/ConformerRigidDockingEngine.h"
#include "Engines/ScoringFunctions/VinaLikeScoringFunction.h"
#include <Engines/Internals/InternalsUtilityFunctions.h>


#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>

#include <Utilities/IntermediateConformerCollector.h>
#include <Structures/InputPostProcessors/VinaCompatibilityPostProcessor.h>

SmolDock::IntermediateConformerCollector* conformerCollector;


int main() {


    /* Setting up the logger */
    boost::log::core::get()->set_filter
            (
#ifdef NDEBUG
            boost::log::trivial::severity >= boost::log::trivial::info
#else
            boost::log::trivial::severity >= boost::log::trivial::debug
#endif
    );

    auto console_logger = boost::log::add_console_log(std::cout);
    console_logger->set_formatter([](boost::log::record_view const &rec, boost::log::formatting_ostream &strm) {
        if (rec[boost::log::trivial::severity] == boost::log::trivial::trace) {
            strm << " T  "; //         use TRACE_LOG(); macro for auto file:line:function
        } else if (rec[boost::log::trivial::severity] == boost::log::trivial::debug) {
            strm << "{D} ";
        } else if (rec[boost::log::trivial::severity] == boost::log::trivial::info) {
            strm << "    ";
        } else if (rec[boost::log::trivial::severity] == boost::log::trivial::warning) {
            strm << "[!] ";
        } else if (rec[boost::log::trivial::severity] >= boost::log::trivial::error) {
            strm << "[E] ";
        }

        strm << rec[boost::log::expressions::smessage];
    });


    BOOST_LOG_TRIVIAL(info) << "SmolDock v0.1";

    std::vector<std::shared_ptr<SmolDock::InputPostProcessor::InputPostProcessor>> postProcessors;
    postProcessors.push_back(std::make_shared<SmolDock::InputPostProcessor::VinaCompatibilityPostProcessor>());

    SmolDock::Protein prot;
    // prot.populateFromPDB("1dpx.pdb"); // Lysozyme
    prot.populateFromPDB("../DockingTests/COX2_Ibuprofen/4PH9_COX2_without_Ibuprofen.pdb", postProcessors); // COX-2
    //prot.populateFromPDB("../DockingTests/Tripeptide/Tripeptide.pdb"); // A random tripeptide

    SmolDock::Molecule mol;
    //mol.populateFromSMILES("CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"); // Ibuprofen
    mol.populateFromPDB("../DockingTests/COX2_Ibuprofen/VINA_Cox2_BestRes.pdb",
                        "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O", /* SMILES hint for bond order*/
                        120 /* seed */,
                        postProcessors);


    SmolDock::iConformer conf_init = mol.getInitialConformer();
    SmolDock::iProtein iprot = prot.getiProtein();
    BOOST_LOG_TRIVIAL(debug) << "Score without docking : "
                             << SmolDock::Score::vina_like_rigid_inter_scoring_func(conf_init,
                                                                                    SmolDock::iTransformIdentityInit(),
                                                                                    iprot);


    SmolDock::PDBWriter pwriter;
    SmolDock::IntermediateConformerCollector collector(&mol, &pwriter);
    conformerCollector = &collector;


    SmolDock::Engine::ConformerRigidDockingEngine docker(10, /* Number of conformer */
                                                         &prot,
                                                         &mol,
                                                         SmolDock::Score::ScoringFunctionType::VinaRigid, /* Scoring function */
                                                         SmolDock::Heuristics::GlobalHeuristicType::RandomRestart, /* Global heuristic */
                                                         SmolDock::Optimizer::LocalOptimizerType::L_BFGS, /* Local optimizer */
                                                         1244 /* Random seed */);

    //docker.setDockingBox(SmolDock::Engine::AbstractDockingEngine::DockingBoxSetting::everything);

    if (!docker.setupDockingEngine()) {
        BOOST_LOG_TRIVIAL(error) << "Error while setting up engine";
        return 2;
    }


    docker.runDockingEngine();


    std::shared_ptr<SmolDock::DockingResult> res = docker.getDockingResult();


    if (res->ligandPoses.empty()) {
        BOOST_LOG_TRIVIAL(error) << "No result to export";
        return 3;
    }


    for (auto &mol: res->ligandPoses) {
        pwriter.addLigand(mol);
    }


    pwriter.writePDB("res.pdb");


    return 0;



    /*
    SmolDock::DockingResultPrinter printer(res);

    printer.printToConsole();
    */
    return 0;

}