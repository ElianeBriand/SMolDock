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

#include <mpi.h>

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
#include <Utilities/AdvancedErrorHandling.h>

#include <thread>

#include <boost/lexical_cast.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>

#include <boost/filesystem.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>


using boost::lexical_cast;

namespace mpi = boost::mpi;
namespace sd = SmolDock;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {

    setupAdvancedErrorHandling();

    tbb::task_scheduler_init tbbInit(std::thread::hardware_concurrency());


    mpi::environment env(argc, argv);
    mpi::communicator world;

    po::options_description desc("Supported options");
    desc.add_options()
            ("help", "help message")
            ("receptor", po::value<std::string>(), "input receptor PDB file")
            ("ligands_csv", po::value<std::string>(), "input ligand CSV. SMILES and Score column mandatory. Comma separated.")
            ("anchor_ligand", po::value<std::string>(), "input anchor ligand mol2 file.")
            ("output_prefix", po::value<std::string>(), "machine readable output file prefix (will save multiple like PREFIX.1.dat)")
            ("log", po::value<std::string>(), "write log to file (in addition to MPI standard output)")
            ("resume", "Attemp resume from last output file (last = based on filename)")
            ("thread_per_node", po::value<int>(), "Thread per MPI worker node")
            ("center_x", po::value<double>(), "Center of search space : x coordinate")
            ("center_y", po::value<double>(), "Center of search space : y coordinate")
            ("center_z", po::value<double>(), "Center of search space : z coordinate")
            ("radius", po::value<double>(), "Search space radius")
            ("batch_size", po::value<int>(), "Optimization algorithm batch size")
            ("num_starting_conformer", po::value<int>(), "Num starting conformers. NB: docking is still flexible, reduce noise")
            ("num_retry", po::value<int>(), "Num retry for each conformer. (Also noice reduction)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    int thread_per_node = -1;
    std::string logfile_path;

    if (vm.count("thread_per_node") != 0) {
        int tpn_ = vm["thread_per_node"].as<int>();
        if(tpn_ > 0) {
            thread_per_node = tpn_;
        } else {
            BOOST_LOG_TRIVIAL(error) << "Number of thread per node must be positive (" << tpn_ << " received)";
            MPI_Abort(MPI_COMM_WORLD, 31);
        }
    }
    
    if (world.rank() == 0) {
        sd::setupLogPrinting(false, false, "[MPI_0] ");

        if (vm.count("log") != 0) {
            logfile_path = vm["log"].as<std::string>();
            sd::logToFile(logfile_path);
            BOOST_LOG_TRIVIAL(info) << "\n\n-- Logging to " << logfile_path <<" --\n\n";
        }
        


        BOOST_LOG_TRIVIAL(info) << "";
        BOOST_LOG_TRIVIAL(info) << "";
        BOOST_LOG_TRIVIAL(info) << "SmolDock v0.2 - MPI Scoring function calibrator";
        BOOST_LOG_TRIVIAL(info) << "Copyright (C) 2019  Eliane Briand";
        BOOST_LOG_TRIVIAL(info) << "    This program comes with ABSOLUTELY NO WARRANTY.";
        BOOST_LOG_TRIVIAL(info) << "    This is free software, and you are welcome to redistribute it.";
        BOOST_LOG_TRIVIAL(info) << "    You should have received a copy of the GNU General Public License version 3";
        BOOST_LOG_TRIVIAL(info) << "    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n\n";

        std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        BOOST_LOG_TRIVIAL(info) << "Starting at " << std::ctime( &now);
        BOOST_LOG_TRIVIAL(info) << "Detected cores on rank 0 node: " << std::thread::hardware_concurrency();
        if(thread_per_node == -1) {
            BOOST_LOG_TRIVIAL(info) << "Worker rank thread policy: one thread per detected core";
        } else {
            BOOST_LOG_TRIVIAL(info) << "Worker rank thread policy: " << thread_per_node << "threads per node";
        }

        std::string receptor_path;
        std::string ligand_csv_path;
        std::string output_prefix = "mpicalibrator_out";
        std::string anchor_ligand_path;
        bool resume = false;
        double center_x, center_y, center_z, radius;
        int batch_size = 10;
        sd::Calibration::MPICalibratorDirector_ResumeObject resumeObject;
        int num_conformer_per_ligand = 5;
        int num_retry_per_conformer = 5;

        if (vm.count("help")) {
            std::cout << desc << "\n";
            MPI_Abort(MPI_COMM_WORLD, 10);
        }

        if (vm.count("receptor") == 0) {
            BOOST_LOG_TRIVIAL(error)
                << "A receptor PDB file path must be provided with option --receptor (see also --help)";
            std::cout << desc << "\n";
            MPI_Abort(MPI_COMM_WORLD, 11);
        } else {
            receptor_path = vm["receptor"].as<std::string>();
            if (!boost::filesystem::exists(receptor_path)) {
                BOOST_LOG_TRIVIAL(error)
                    << "Receptor file " << receptor_path << " does not exist.";
                MPI_Abort(MPI_COMM_WORLD, 21);
            }
        }

        if (vm.count("ligands_csv") == 0) {
            BOOST_LOG_TRIVIAL(error) << "A ligand CSV must be provided with option --ligands_csv (see also --help). SMILES and Score column mandatory. Comma separated.";
            MPI_Abort(MPI_COMM_WORLD, 12);
        } else {
            ligand_csv_path = vm["ligands_csv"].as<std::string>();
            if (!boost::filesystem::exists(ligand_csv_path)) {
                BOOST_LOG_TRIVIAL(error)
                    << "Ligand CSV file " << ligand_csv_path << " does not exist.";
                MPI_Abort(MPI_COMM_WORLD, 22);
            }
        }

        if (vm.count("anchor_ligand") == 0) {
            BOOST_LOG_TRIVIAL(error) << "An anchor ligand mol2 file must be provided with option --ligands_csv (see also --help). SMILES and Score column mandatory. Comma separated.";
            MPI_Abort(MPI_COMM_WORLD, 16);
        } else {
            anchor_ligand_path = vm["anchor_ligand"].as<std::string>();
            if (!boost::filesystem::exists(anchor_ligand_path)) {
                BOOST_LOG_TRIVIAL(error)
                    << "Anchor ligand mol2 file " << anchor_ligand_path << " does not exist.";
                MPI_Abort(MPI_COMM_WORLD, 26);
            }
        }


        if (vm.count("output_prefix") == 0) {
            BOOST_LOG_TRIVIAL(error) << "An output prefix must be provided with option --output_prefix (see also --help). ";
            MPI_Abort(MPI_COMM_WORLD, 13);
        } else {
            output_prefix = vm["output_prefix"].as<std::string>();
            BOOST_LOG_TRIVIAL(info)
                    << "Using output prefix " << output_prefix << " for this run.";
        }

        if ((vm.count("center_x") == 0) ||
            (vm.count("center_y") == 0) ||
            (vm.count("center_z") == 0) ||
            (vm.count("radius") == 0) ) {
            BOOST_LOG_TRIVIAL(error) << "Parameters for the search volume must be defined using --center_x,--center_y,--center_z,--radius (see also --help). ";
            MPI_Abort(MPI_COMM_WORLD, 14);
        } else {
            center_x = vm["center_x"].as<double>();
            center_y = vm["center_y"].as<double>();
            center_z = vm["center_z"].as<double>();
            radius = vm["radius"].as<double>();

            BOOST_LOG_TRIVIAL(info)
                << "Using search volume parameters: center [" << center_x << ", " << center_y << ", " << center_z <<
                            " ] r = " << radius << " A for this run.";
        }

        if (vm.count("resume") != 0) {
            BOOST_LOG_TRIVIAL(info) << "Attempting resume...";
            const boost::regex my_filter(output_prefix + std::string("([0-9]+)\\.restore-state") );

            int max_num = -1;
            boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
            for( boost::filesystem::directory_iterator i("./"); i != end_itr; ++i )
            {
                // Skip if not a file
                if( !boost::filesystem::is_regular_file( i->status() ) ) continue;

                boost::smatch what;

                // Skip if no match for V2:
                //if( !boost::regex_match( i->leaf(), what, my_filter ) ) continue;
                // For V3:
                std::string path_as_string = i->path().filename().string();
                bool res = boost::regex_match( path_as_string, what, my_filter, boost::match_extra) ;
                if(!res) continue;
                int numFile = lexical_cast<int>(what[1]);
                if(numFile > max_num) {
                    max_num = numFile;
                }
            }

            if(max_num != -1) {
                std::string filename = output_prefix + std::to_string(max_num) + std::string(".restore-state");
                BOOST_LOG_TRIVIAL(info) << "Loading from: " << filename;
                std::ifstream ifs(filename);
                boost::archive::text_iarchive ia(ifs);
                ia >> resumeObject;
                resume = true;
                BOOST_LOG_TRIVIAL(info) << "State loaded successfully";
            }else {
                BOOST_LOG_TRIVIAL(info) << "Resume cannot proceed without a restore-state file";
                MPI_Abort(MPI_COMM_WORLD, 41);
            }


        }


        if (vm.count("batch_size") != 0) {
            unsigned int batch_size_ = vm["batch_size"].as<int>();
            if(batch_size_ <= 0) {
                BOOST_LOG_TRIVIAL(info) << "Batch size must be positive integer";
                MPI_Abort(MPI_COMM_WORLD, 26);
            }
            batch_size = batch_size_;
        }


        if (vm.count("num_starting_conformer") != 0) {
            unsigned int num_conformer_per_ligand_ = vm["num_starting_conformer"].as<int>();
            if(num_conformer_per_ligand_ <= 0) {
                BOOST_LOG_TRIVIAL(info) << "Number of starting conformer must be positive";
                MPI_Abort(MPI_COMM_WORLD, 27);
            }
            num_conformer_per_ligand = num_conformer_per_ligand_;
        }

        if (vm.count("num_retry") != 0) {
            unsigned int num_retry_per_conformer_ = vm["num_retry"].as<int>();
            if(num_retry_per_conformer_ <= 0) {
                BOOST_LOG_TRIVIAL(info) << "Number of retry per conformer must be positive";
                MPI_Abort(MPI_COMM_WORLD, 28);
            }
            num_retry_per_conformer = num_retry_per_conformer_;
        }



        auto cdirector = std::make_shared<sd::Calibration::MPICalibratorDirector>(env, world,
                                                                                  sd::Score::ScoringFunctionType::VinaCovalentReversible,
                                                                                  sd::Heuristics::GlobalHeuristicType::SimulatedAnnealing,
                                                                                  sd::Optimizer::LocalOptimizerType::L_BFGS,
                                                                                  1000, // Max epoch
                                                                                  0.5e-7, // Initial learning rate
                                                                                  32574, // RNG seed
                                                                                  num_conformer_per_ligand, // num generated starting conformer per ligand
                                                                                  num_retry_per_conformer, // num retry per conformer (best score is kept)
                                                                                  batch_size, //batch size
                                                                                  sd::Heuristics::emptyParameters);

        if(resume) {
            BOOST_LOG_TRIVIAL(debug) << "Attempting state restore";
            cdirector->restoreResumeState(resumeObject);
            BOOST_LOG_TRIVIAL(info) << "State restored successfully";
        }else {
            cdirector->coefficientsToCalibrate({"Gauss1","Gauss2","RepulsionExceptCovalent","Hydrophobic","Hydrogen"});
            //cdirector->coefficientsToCalibrate({"CovalentReversible"});
        }




        sd::Engine::AbstractDockingEngine::DockingBoxSetting setting;
        setting.type = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Type::centeredAround;
        setting.center = {center_x, center_y, center_z};
        setting.radius = radius;
        sd::Calibration::Calibrator::ReceptorID recID1 =
                cdirector->addReceptorFromFile(receptor_path, setting);


        sd::CSVReader chembl_csv(ligand_csv_path,",",true, "µ");


        std::vector<std::map<std::string,std::string>> chembl_data = chembl_csv.getRowsAsMap();


        for (unsigned int j = 0; j < chembl_data.size(); ++j) {
            std::map<std::string,std::string> row = chembl_data.at(j);
            std::string smiles_str = row["SMILES"];
            std::string score_str = row["Score"];
            double score = boost::lexical_cast<double>(score_str);
            cdirector->addReferenceLigand_SMILES_deltaG(recID1,smiles_str, score);
        }


        cdirector->setupCalibration();
        cdirector->runCalibration();

    } else {
        if(thread_per_node == -1) {
            thread_per_node = std::thread::hardware_concurrency();
        }
        tbb::task_scheduler_init tbbInit(thread_per_node);
        //tbb::task_scheduler_init tbbInit(1);
        sd::setupLogPrinting(false, false, "[MPI_" + lexical_cast<std::string>(world.rank()) +"] ");
        auto cnode = std::make_shared<sd::Calibration::MPICalibratorNode>(env, world);
        cnode->runNode();
    }


    return 0;


}
