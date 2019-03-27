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
#include <Frontends/FrontendCommon.h>

#include <boost/log/trivial.hpp>

#define USE_BOOST_KARMA
#include <bprinter/table_printer.h>

#include <tbb/task_scheduler_init.h>

namespace sd = SmolDock;



int main() {

    //tbb::task_scheduler_init tbbInit(2);

    /* Setting up the logger */
    sd::setupLogPrinting();


    BOOST_LOG_TRIVIAL(info) << "SmolDock " << sd::getVersionString();

    std::vector<std::shared_ptr<sd::InputModifier::InputModifier>> modifiers;
    modifiers.push_back(std::make_shared<sd::InputModifier::VinaCompatibility>());

    bool succeeded = true;

    sd::Protein prot;
    succeeded &= prot.populateFromPDB("../DockingTests/COX2_Ibuprofen/4PH9_COX2_without_Ibuprofen.pdb", modifiers); // COX-2


    sd::Molecule mol;
    succeeded &= mol.populateFromMol2File("../DockingTests/COX2_Ibuprofen/VINA_Cox2_BestRes_Charged.mol2", 120, modifiers); // Ibuprofen


    if (!succeeded) {
        BOOST_LOG_TRIVIAL(error) << "Error while loading input files";
        return 2;
    }



    sd::iConformer conf_init = mol.getInitialConformer();
    sd::iProtein iprot = prot.getiProtein();
    sd::iTransform tr = sd::iTransformIdentityInit();

    double scoreWithoutDocking = sd::Score::vina_like_rigid_inter_scoring_func(conf_init,
                                                                                     tr,
                                                                                     iprot);
    BOOST_LOG_TRIVIAL(debug) << "Score without docking : "
                             << scoreWithoutDocking;


    sd::PDBWriter pwriter;


    sd::Engine::ConformerDockingEngine docker(5, 4,&prot,&mol,
                                                       sd::Score::ScoringFunctionType::Vina,
                                                       sd::Heuristics::GlobalHeuristicType::SimulatedAnnealing,
                                                       sd::Optimizer::LocalOptimizerType::L_BFGS,
                                                       1244);


    sd::Engine::AbstractDockingEngine::DockingBoxSetting setting;
    setting.type = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround;
    setting.center = {10.0, 22.0, 25.0};
    setting.radius = 10.0;

    docker.setDockingBox(setting);

    if (!docker.setupDockingEngine()) {
        BOOST_LOG_TRIVIAL(error) << "Error while setting up engine";
        return 2;
    }



    docker.runDockingEngine();

    BOOST_LOG_TRIVIAL(debug) << "Reminder : Score without docking : "
                             << scoreWithoutDocking;

    std::shared_ptr<sd::DockingResult> res = docker.getDockingResult();

    if (res->ligandPoses.empty()) {
        BOOST_LOG_TRIVIAL(error) << "No result to export";
        return 3;
    }


    for (auto &mol1: res->ligandPoses) {
        pwriter.addLigand(mol1);
    }


    pwriter.writePDB("res_speedtest.pdb");

    return 0;

////
////    if (!succeeded) {
////        BOOST_LOG_TRIVIAL(error) << "Error while loading input files";
////        return 2;
////    }
////
////    sd::iConformer conf_init = mol.getInitialConformer();
////    sd::iProtein iprot = prot.getiProtein();
////
////    double scoreWithoutDocking = sd::Score::vina_like_rigid_inter_scoring_func(conf_init,
////                                                                                     sd::iTransformIdentityInit(),
////                                                                                     iprot);
////    BOOST_LOG_TRIVIAL(debug) << "Score without docking : "
////                             << scoreWithoutDocking;
////
////
////    sd::PDBWriter pwriter;
////
////
/////*
////     *         struct Parameters {
////            unsigned int maxIterations = 10000;
////            double initTemp = 10000.0;
////            unsigned int initialNoTempDropMoves = 1000;
////            unsigned int moveCtrlSweep = 100;
////            double tolerance = 1e-3;
////            unsigned int maxToleranceSweep = 3;
////            double maxMoveSize = 5.0;
////            double initMoveSize = 0.3;
////            double gain = 0.3;
////        };
////     * */
////
////
//////    sd::Heuristics::SimulatedAnnealing::Parameters param;
////    sd::Heuristics::HeuristicParameters param((std::in_place_type_t<sd::Heuristics::SimulatedAnnealing::Parameters>()));
////
////    std::vector<std::tuple<double,double,double,double,double,double>> gridSearchRes;
////
////    for(unsigned int iteration = 1; iteration < 20; iteration++)
////    {
////
////        double paramValue = iteration;
////        std::get<sd::Heuristics::SimulatedAnnealing::Parameters>(param).maxMoveSize = paramValue;
////
////        sd::Engine::ConformerDockingEngine docker(5,
////                                                       3,
////                                                       &prot,
////                                                       &mol,
////                                                       sd::Score::ScoringFunctionType::Vina,
////                                                       sd::Heuristics::GlobalHeuristicType::SimulatedAnnealing,
////                                                       sd::Optimizer::LocalOptimizerType::L_BFGS,
////                                                       1244,
////                                                       param);
////
////
////        sd::Engine::AbstractDockingEngine::DockingBoxSetting setting;
////        setting.type = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround;
////        setting.center = {10.0, 22.0, 25.0};
////        setting.radius = 10.0;
////
////        docker.setDockingBox(setting);
////
////        if (!docker.setupDockingEngine()) {
////            BOOST_LOG_TRIVIAL(error) << "Error while setting up engine";
////            return 2;
////        }
////
////
////        docker.runDockingEngine();
////
////        BOOST_LOG_TRIVIAL(debug) << "Reminder : Score without docking : "
////                                 << scoreWithoutDocking;
////
////        double meanDuration, stdDevDuration, meanScore, stdDevScore, bestScore;
////        std::tie(meanDuration,stdDevDuration) = docker.getMeanStdDevDuration();
////        std::tie(meanScore,stdDevScore) = docker.getMeanStdDevScore();
////        bestScore = docker.getBestScore();
////
////        gridSearchRes.emplace_back(std::make_tuple(paramValue,
////                                                meanDuration, stdDevDuration,bestScore, meanScore, stdDevScore));
////    }
////
////    bprinter::TablePrinter tp(&std::cout);
////    tp.AddColumn("maxMoveSize",22);
////    tp.AddColumn("Mean Duration", 22);
////    tp.AddColumn("StdDev Duration", 22);
////    tp.AddColumn("Best Score", 22);
////    tp.AddColumn("Mean Score", 22);
////    tp.AddColumn("StdDev Score", 22);
////
////    tp.PrintHeader();
////    for(auto& record: gridSearchRes)
////    {
////        tp <<  std::get<0>(record)
////            << std::get<1>(record)
////            << std::get<2>(record)
////            << std::get<3>(record)
////            << std::get<4>(record)
////            << std::get<5>(record);
////    }
////    tp.PrintFooter();
//
//
//
//
///*
//    std::shared_ptr<sd::DockingResult> res = docker.getDockingResult();
//
//
//    if (res->ligandPoses.empty()) {
//        BOOST_LOG_TRIVIAL(error) << "No result to export";
//        return 3;
//    }
//
//
//    for (auto &mol1: res->ligandPoses) {
//        pwriter.addLigand(mol1);
//    }
//
//
//    pwriter.writePDB("res.pdb");
//
//    */
//
//
//    return 0;
//
//
//
//    /*
//    sd::DockingResultPrinter printer(res);
//
//    printer.printToConsole();
//    */
//    return 0;

}
