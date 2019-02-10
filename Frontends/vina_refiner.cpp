/*
 * Copyright (c) 2019 Eliane Briand
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
#include <string>

#include "Structures/Molecule.h"
#include "Structures/Protein.h"
#include "Engines/PoseRefiner.h"

#include "Structures/InputPostProcessors/VinaCompatibilityPostProcessor.h"

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <Utilities/PDBWriter.h>


namespace po = boost::program_options;


int main(int argc, char* argv[]) {

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


    BOOST_LOG_TRIVIAL(info) << "vina_refiner utility (SmolDock v0.1)";
    BOOST_LOG_TRIVIAL(info) << "Copyright (C) 2019  Eliane Briand";
    BOOST_LOG_TRIVIAL(info) << "    This program comes with ABSOLUTELY NO WARRANTY.";
    BOOST_LOG_TRIVIAL(info) << "    This is free software, and you are welcome to redistribute it.";
    BOOST_LOG_TRIVIAL(info) << "    You should have received a copy of the GNU General Public License version 3";
    BOOST_LOG_TRIVIAL(info) << "    along with this program.  If not, see <https://www.gnu.org/licenses/>.\n";


    po::options_description desc("Supported options");
    desc.add_options()
            ("help", "help message")
            ("receptor", po::value< std::string >(), "input receptor PDB file")
            ("ligand", po::value< std::string >(), "input ligand PDB file")
            ("output", po::value< std::string >(), "output PDB file for ligand conformation")
            ("log", po::value< std::string >(), "write log to file (in addition to standard output)")
            ("force", "overwrite output file if it already exist")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        BOOST_LOG_TRIVIAL(info) << desc << "\n";
        return 1;
    }

    bool argument_error = false;
    if (vm.count("receptor") == 0)
    {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error) << "A receptor PDB file path must be provided with option --receptor (see also --help)";
             //<< vm["include-path"].as< std::string> >() << "\n";
    }

    if (vm.count("ligand") == 0)
    {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error) << "A ligand PDB file path must be provided with option --ligand (see also --help)";
             //<< vm["input-file"].as< vector<string> >() << "\n";
    }

    if (vm.count("output") == 0)
    {
        argument_error = true;
        BOOST_LOG_TRIVIAL(error) << "An output ligand PDB file path must be provided with option --output (see also --help)";
        //<< vm["input-file"].as< vector<string> >() << "\n";
    }

    if(argument_error == true)
        return 1;

    std::string path_to_input_receptor = vm["receptor"].as< std::string >();
    std::string path_to_input_ligand = vm["ligand"].as< std::string >();
    std::string path_to_output_ligand = vm["output"].as< std::string >();

    bool overwriteFiles = false;
    if(vm.count("force"))
        overwriteFiles = true;

    bool produceLog = false;
    std::string path_to_log;
    if (vm.count("log"))
    {
        produceLog = true;
        path_to_log = vm["log"].as< std::string>();
        auto file_logger = boost::log::add_file_log
                (
                        boost::log::keywords::file_name = path_to_log
                );

        file_logger->set_formatter([](boost::log::record_view const &rec, boost::log::formatting_ostream &strm) {
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

    }



    BOOST_LOG_TRIVIAL(info) << "input receptor path : " << path_to_input_receptor;
    BOOST_LOG_TRIVIAL(info) << "input ligand path   : " << path_to_input_ligand;
    BOOST_LOG_TRIVIAL(info) << "output ligand path  : " << path_to_output_ligand;
    if(produceLog)
        BOOST_LOG_TRIVIAL(info) << "log file path       : " << path_to_log;

    if ( !boost::filesystem::exists( path_to_input_receptor ) )
    {
        BOOST_LOG_TRIVIAL(error) << "Receptor file not found : " << path_to_input_receptor;
        return 1;
    }


    if ( !boost::filesystem::exists( path_to_input_ligand ) )
    {
        BOOST_LOG_TRIVIAL(error) << "Ligand file not found : " << path_to_input_receptor;
        return 1;
    }


    if ( boost::filesystem::exists( path_to_output_ligand ) )
    {
        if(overwriteFiles)
            BOOST_LOG_TRIVIAL(info) << "Output ligand file already exists, will be overwritten.";
        else
        {
            BOOST_LOG_TRIVIAL(info) << "Output file already exist, aborting";
            return 1;
        }
    }




    std::vector<std::shared_ptr<SmolDock::InputPostProcessor::InputPostProcessor>> postProcessors;
    postProcessors.push_back(std::make_shared<SmolDock::InputPostProcessor::VinaCompatibilityPostProcessor>());

    SmolDock::Protein prot;
    prot.populateFromPDB(path_to_input_receptor, postProcessors); // COX-2

    SmolDock::Molecule mol;
    mol.populateFromPDB(path_to_input_ligand,"", /* No SMILES hint for bond order*/
                        120 /* seed */,
                        postProcessors);



    SmolDock::Engine::PoseRefiner vinaPoseRefiner(&prot,&mol,
                                          SmolDock::Score::ScoringFunctionType::VinaRigid,
                                          SmolDock::Optimizer::LocalOptimizerType::L_BFGS,120);

    vinaPoseRefiner.refinePose();
    vinaPoseRefiner.applyToLigand();

    double init_score = vinaPoseRefiner.getInitialScore();
    double final_score = vinaPoseRefiner.getFinalScore();
    double score_diff = vinaPoseRefiner.getScoreDifference();


    SmolDock::PDBWriter w;
    w.addLigand(mol);

    BOOST_LOG_TRIVIAL(info) << "Initial score : " << std::fixed << std::setprecision(6)<< init_score;
    BOOST_LOG_TRIVIAL(info) << "Final score   : " << std::fixed << std::setprecision(6)<< final_score;
    BOOST_LOG_TRIVIAL(info) << "Improvement   :  " << std::fixed << std::setprecision(6)<< score_diff;

    w.writePDB(path_to_output_ligand);

    return 0;
}