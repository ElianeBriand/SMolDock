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
#include <optional>

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



std::optional<sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape> boxShapeFromString(const std::string& s);

int main(int argc, char *argv[]) {

    setupAdvancedErrorHandling();
    sd::setupLogPrinting();

    BOOST_LOG_TRIVIAL(info) << ""; 
    BOOST_LOG_TRIVIAL(info) << "SmolDock v0.2";
    BOOST_LOG_TRIVIAL(info) << "Copyright (C) 2019  Eliane Briand";
    BOOST_LOG_TRIVIAL(info) << "    This program comes with ABSOLUTELY NO WARRANTY.";
    BOOST_LOG_TRIVIAL(info) << "    This is free software, and you are welcome to redistribute it.";
    BOOST_LOG_TRIVIAL(info) << "    You should have received a copy of the GNU General Public License version 3";
    BOOST_LOG_TRIVIAL(info) << "    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n\n";


    po::options_description general_opt("General options");
    general_opt.add_options()
            ("help", "help message");

    po::options_description input_opt("Input options");
    input_opt.add_options()
            ("receptor", po::value<std::string>(), "input receptor PDB file")
            ("ligand", po::value<std::string>(), "input ligand mol or mol2 file")
            ;

    po::options_description output_opt("Output options");
    output_opt.add_options()
            ("output", po::value<std::string>(), "output PDB file for ligand conformation")
            ("log", po::value<std::string>(), "write log to file (in addition to standard output)")
            ("force", "overwrite output file if it already exist")
            ;
    
    po::options_description searchvol_opt("Search volume options");
    searchvol_opt.add_options()
            ("center_x", po::value<double>(), "Center of search volume : x coordinate")
            ("center_y", po::value<double>(), "Center of search volume : y coordinate")
            ("center_z", po::value<double>(), "Center of search volume : z coordinate")
            ("boxtype", po::value<std::string>(), "Type of search volume [sphere, box, whole_protein]")
            ("radius", po::value<double>(), "Search volume radius (if spherical)")
            ("size_x", po::value<double>(), "Search space box side: x axis")
            ("size_y", po::value<double>(), "Search space box side: y axis")
            ("size_z", po::value<double>(), "Search space box side: z axis")
            ;

    po::options_description perf_opt("Performance and exhaustiveness options");
    perf_opt.add_options()
            ("cpu", po::value<int>(), "Number of thread")
            ("exhaustiveness", po::value<int>(), "Holistic search exhaustiveness parameter (higher is slower but more diligent)")
            ;

    po::options_description misc_opt("Misc options");
    misc_opt.add_options()
            ("seed", po::value<unsigned int>()->default_value(125), "RNG seed")
            ;


    po::options_description cmdline_options;
    cmdline_options.add(general_opt)
            .add(input_opt)
            .add(output_opt)
            .add(searchvol_opt)
            .add(perf_opt)
            .add(misc_opt);

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
    po::notify(vm);

    if (vm.count("help")) {
        BOOST_LOG_TRIVIAL(info) << cmdline_options << "\n";
        return 1;
    }

    bool argument_error = false;
    if (vm.count("receptor") == 0) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error)
            << "A receptor PDB file path must be provided with option --receptor";
    }

    if (vm.count("ligand") == 0) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error) << "A ligand mol or mol2 file path must be provided with option --ligand";
    }

    if (vm.count("output") == 0) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error)
            << "An output ligand PDB file path must be provided with option --output";
    }

    if ((vm.count("center_x") == 0) ||
        (vm.count("center_y") == 0) ||
        (vm.count("center_z") == 0) ||
        (vm.count("boxtype") == 0)  || 
        ( (vm.count("radius") == 0) && 
            ((vm.count("size_x") == 0) || (vm.count("size_y") == 0) || (vm.count("size_z") == 0))
         )
     ) {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error)
            << "Search volume parameters must be specified:";
        BOOST_LOG_TRIVIAL(error)
            << "    center + radius + sphere, or center + size + box";
    }

    if (argument_error == true) {
        BOOST_LOG_TRIVIAL(info) << ""; 
        BOOST_LOG_TRIVIAL(info) << cmdline_options << "\n"; 
        return 1;
    }

    unsigned int seed = vm["seed"].as<unsigned int>();

    std::string path_to_input_receptor = vm["receptor"].as<std::string>();
    std::string path_to_input_ligand = vm["ligand"].as<std::string>();
    std::string path_to_output_ligand = vm["output"].as<std::string>();

    double center_x = vm["center_x"].as<double>();
    double center_y = vm["center_y"].as<double>();
    double center_z = vm["center_z"].as<double>();

    sd::Engine::AbstractDockingEngine::DockingBoxSetting setting;
    auto shape_opt = boxShapeFromString(vm["boxtype"].as<std::string>());
    if(!shape_opt.has_value())
    {
        BOOST_LOG_TRIVIAL(error)
            << "Unknown box type: " << vm["boxtype"].as<std::string>();
        return 0;
    }
    setting.shape = shape_opt.value();
    setting.center = {center_x, center_y, center_z};
    if(setting.shape == sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::cube) {
        double size_x = vm["size_x"].as<double>();
        double size_y = vm["size_y"].as<double>();
        double size_z = vm["size_z"].as<double>();
        setting.size = {size_x, size_y, size_z};
    }else if (setting.shape == sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::sphere) {
        double r = vm["radius"].as<double>();
        setting.radius = r;
    }
    


    bool overwriteFiles = false;
    if (vm.count("force"))
        overwriteFiles = true;

    bool produceLog = false;
    std::string path_to_log;
    if (vm.count("log")) {
        produceLog = true;
        path_to_log = vm["log"].as<std::string>();
        sd::logToFile(path_to_log);
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
                                seed,
                                postProcessors);
    }else if (pathLigand.extension() == ".mol2") {
        mol.populateFromMol2File(pathLigand.string(),
                                 seed,
                                 postProcessors);
    }else {
        BOOST_LOG_TRIVIAL(error) << "Unsupported ligand file extension \"" << pathLigand.extension() << "\"";
        BOOST_LOG_TRIVIAL(info) << "Supported extension : mol mol2";
    }


    SmolDock::Engine::ConformerDockingEngine docker(20, /* Number of conformer */
                                                    10, /* Retry per conformer */
                                                    &prot,
                                                    &mol,
                                                    SmolDock::Score::ScoringFunctionType::Vina, /* Scoring function */
                                                    SmolDock::Heuristics::GlobalHeuristicType::SimulatedAnnealing, /* Global heuristic */
                                                    SmolDock::Optimizer::LocalOptimizerType::L_BFGS, /* Local optimizer */
                                                    seed /* Random seed */);

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


std::optional<sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape> boxShapeFromString(const std::string& s) {
    std::optional<sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape> ret;

    if(s == "cube") {
        ret = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::cube;
    }else if(s == "sphere"){
        ret = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::sphere;
    }else if(s == "whole_protein"){
        ret = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::whole_protein;
    }

    return ret;

}