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
#include <chrono>
#include <ctime>

#include <boost/log/trivial.hpp>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


#include <Structures/Molecule.h>
#include <Structures/Protein.h>
#include <Structures/Atom.h>
#include <Structures/InputModifiers/InputModifierInterface.h>
#include <Structures/InputModifiers/VinaCompatibility.h>

#include <Engines/ConformerDockingEngine.h>
#include <Engines/ScoringFunctions/VinaLike.h>
#include <Engines/ScoringFunctions/VinaLikeRigid.h>
#include <Engines/ScoringFunctions/VinaLikeCommon.h>
#include <Engines/Internals/InternalsUtilityFunctions.h>

#include <Utilities/DockingResultPrinter.h>
#include <Utilities/PDBWriter.h>
#include <Utilities/CSVReader.h>
#include <Utilities/SMARTSMatcher.h>
#include <Utilities/Version.h>
#include <Utilities/AdvancedErrorHandling.h>

#include <Utilities/Calibration/Calibrator.h>

#include <Frontends/FrontendCommon.h>


#define USE_BOOST_KARMA
#include <bprinter/table_printer.h>

#include <tbb/task_scheduler_init.h>



namespace sd = SmolDock;
namespace po = boost::program_options;



int main(int argc, char *argv[]) {

    setupAdvancedErrorHandling();
    sd::setupLogPrinting();


    BOOST_LOG_TRIVIAL(info) << "\n\nSmolDock v0.2";
    BOOST_LOG_TRIVIAL(info) << "Copyright (C) 2019  Eliane Briand";
    BOOST_LOG_TRIVIAL(info) << "    This program comes with ABSOLUTELY NO WARRANTY.";
    BOOST_LOG_TRIVIAL(info) << "    This is free software, and you are welcome to redistribute it.";
    BOOST_LOG_TRIVIAL(info) << "    You should have received a copy of the GNU General Public License version 3";
    BOOST_LOG_TRIVIAL(info) << "    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n\n";


    po::options_description desc("Supported options");
    desc.add_options()
            ("help", "help message")
            ("receptor", po::value<std::string>(), "input receptor PDB file")
            ("ligand", po::value<std::string>(), "input ligand mol or mol2 file")
            ("output", po::value<std::string>(), "output PDB file for ligand conformation")
            ("log", po::value<std::string>(), "write log to file (in addition to standard output)")
            ("force", "overwrite output file if it already exist");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        BOOST_LOG_TRIVIAL(info) << desc << "\n";
        return 1;
    }

    bool argument_error = false;
    if (vm.count("receptor") == 0) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error)
            << "A receptor PDB file path must be provided with option --receptor (see also --help)";
        //<< vm["include-path"].as< std::string> >() << "\n";
    }

    if (vm.count("ligand") == 0) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error) << "A ligand mol or mol2 file path must be provided with option --ligand (see also --help)";
        //<< vm["input-file"].as< vector<string> >() << "\n";
    }

    if (vm.count("output") == 0) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error)
            << "An output ligand PDB file path must be provided with option --output (see also --help)";
        //<< vm["input-file"].as< vector<string> >() << "\n";
    }

    if (argument_error == true)
        return 1;

    std::string path_to_input_receptor = vm["receptor"].as<std::string>();
    std::string path_to_input_ligand = vm["ligand"].as<std::string>();
    std::string path_to_output_ligand = vm["output"].as<std::string>();

    bool overwriteFiles = false;
    if (vm.count("force"))
        overwriteFiles = true;

    bool produceLog = false;
    std::string path_to_log;
    if (vm.count("log")) {
        produceLog = true;
        path_to_log = vm["log"].as<std::string>();
        logToFile(path_to_log);
    }


    BOOST_LOG_TRIVIAL(info) << "input receptor path : " << path_to_input_receptor;
    BOOST_LOG_TRIVIAL(info) << "input ligand path   : " << path_to_input_ligand;
    BOOST_LOG_TRIVIAL(info) << "output ligand path  : " << path_to_output_ligand;
    if (produceLog)
        BOOST_LOG_TRIVIAL(info) << "log file path       : " << path_to_log;

    if (!boost::filesystem::exists(path_to_input_receptor) || !boost::filesystem::is_regular_file(path_to_input_receptor)) {
        BOOST_LOG_TRIVIAL(error) << "Receptor file not found (or not a file) : " << path_to_input_receptor;
        return 1;
    }


    if (!boost::filesystem::exists(path_to_input_ligand) || !boost::filesystem::is_regular_file(path_to_input_ligand)) {
        BOOST_LOG_TRIVIAL(error) << "Ligand file not found (or not a file) : " << path_to_input_receptor;
        return 1;
    }


    if (boost::filesystem::exists(path_to_output_ligand)) {
        if (overwriteFiles)
            BOOST_LOG_TRIVIAL(info) << "Output ligand file already exists, will be overwritten.";
        else {
            BOOST_LOG_TRIVIAL(info) << "Output file already exist, aborting";
            return 1;
        }
    }

    boost::filesystem::path pathOutput(path_to_output_ligand);
    boost::filesystem::path pathReceptor(path_to_input_receptor);
    boost::filesystem::path pathLigand(path_to_input_ligand);

    std::vector<std::shared_ptr<SmolDock::InputModifier::InputModifier>> postProcessors;
    postProcessors.push_back(std::make_shared<SmolDock::InputModifier::VinaCompatibility>());

    SmolDock::Protein prot;
    prot.populateFromPDB(path_to_input_receptor, postProcessors); // COX-2

    SmolDock::Molecule mol;

    if(pathLigand.extension() == ".mol")
    {
        mol.populateFromMolFile(pathLigand.string(),
                                120 /* seed */,
                                postProcessors);
    }else if (pathLigand.extension() == ".mol2") {
        mol.populateFromMol2File(pathLigand.string(),
                                 120 /* seed */,
                                 postProcessors);
    }else {
        BOOST_LOG_TRIVIAL(error) << "Unsupported ligand file extension \"" << pathLigand.extension() << "\"";
        BOOST_LOG_TRIVIAL(info) << "Supported extension : mol mol2";
    }


    SmolDock::Engine::ConformerDockingEngine docker(20, /* Number of conformer */
                                                    10, /* Retry per conformer */
                                                    &prot,
                                                    &mol,
                                                    SmolDock::Score::ScoringFunctionType::VinaRigid, /* Scoring function */
                                                    SmolDock::Heuristics::GlobalHeuristicType::IteratedLocalSearch, /* Global heuristic */
                                                    SmolDock::Optimizer::LocalOptimizerType::L_BFGS, /* Local optimizer */
                                                    1244 /* Random seed */);


    SmolDock::Engine::AbstractDockingEngine::DockingBoxSetting setting;
    setting.type = SmolDock::Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround;
    setting.center = {10.0, 22.0, 25.0};
    setting.radius = 10.0;

    docker.setDockingBox(setting);

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

    SmolDock::PDBWriter pwriter;
    for (auto &mol1: res->ligandPoses) {
        pwriter.addLigand(mol1);
    }


    pwriter.writePDB(path_to_output_ligand);


    BOOST_LOG_TRIVIAL(info) << "Done.";


    return 0;
}
