//
// Created by eliane on 28/11/18.
//

#include <memory>
#include <chrono>
#include <thread>

#include "ConformerRigidDockingEngine.h"
#include "Internals/iConformer.h"

#include "Utilities/TimingsLog.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>


namespace SmolDock {
    namespace Engine {


        ConformerRigidDockingEngine::ConformerRigidDockingEngine(unsigned int conformer_num) {
            this->conformer_num = conformer_num;
            assert(conformer_num != 0);
        }

        bool ConformerRigidDockingEngine::setProtein(Protein *p) {
            this->orig_protein = p;
            return true;
        }

        bool ConformerRigidDockingEngine::setLigand(Molecule *m) {
            this->orig_ligand = m;

            // Getting a handle on the RWMol
            this->rwmol = this->orig_ligand->getInternalRWMol();


            return true;
        }

        bool ConformerRigidDockingEngine::setupDockingEngine() {

            record_timings(begin_setup);

            std::uniform_int_distribution<> dis(0, std::numeric_limits<int>::max());

            record_timings(start_conformersgen);

            this->orig_ligand->generateConformers(this->viConformers, this->conformer_num, dis(this->rnd_generator));

            record_timings(done_conformersgen);




            record_timings(done_setup);


#ifdef SMOLDOCK_VERBOSE_DEBUG
            auto total_setup_time = static_cast< std::chrono::duration<double> >(done_setup - begin_setup).count();
            auto timings_conformer = static_cast< std::chrono::duration<double> >(done_conformersgen - start_conformersgen).count();

            std::cout << "[D] Timings ConformerRigidDockingEngine::setupDockingEngine : "
                                  << total_setup_time
                                  << " seconds total\n      Conformer: "
                                  << timings_conformer
                                  <<"s [n=" << this->conformer_num << "->" << this->viConformers.size() << "] "
                                  << timings_conformer/this->conformer_num << "s each" << std::endl
                                  << "      Other: 0"
                                  <<"s"
                                  << std::endl;
#endif


            return false;
        }

        void ConformerRigidDockingEngine::runDockingEngine() {
            std::cout << "    |";
            for(int i = 0; i < 98; i++)
                std::cout << " ";
            std::cout << "| 100%" << std::endl << "    " ;
            std::cout.flush();
            unsigned int quorum = (this->viConformers.size() - (this->viConformers.size() % 100)) / 100;
            unsigned int current = 0;
            for(auto& conformer : this->viConformers)
            {
                using namespace std::chrono_literals;
                std::this_thread::sleep_for(0.3s);
                if(current % quorum == 0) {
                    std::cout << ".";
                    std::cout.flush();
                }
                current++;
            }
            std::cout << std::endl;

        }

        std::shared_ptr<DockingResult> ConformerRigidDockingEngine::getDockingResult() {
            return std::make_shared<DockingResult>();
        }


        bool ConformerRigidDockingEngine::setDockingBox(AbstractDockingEngine::DockingBoxSetting setting) {
            if (setting != DockingBoxSetting::everything) {
                std::cout << "[!] DockingBoxSetting (that is not DockingBoxSetting::everything) is not yet implemented."
                          << std::endl;
                std::cout << "[ ] Running as if DockingBoxSetting::everything was passed" << std::endl;
                return false;
            }
            return true;
        }


        void ConformerRigidDockingEngine::setRandomSeed(int seed) {
            this->random_seed = seed;
            this->rnd_generator.seed(this->random_seed);
        }


    }

}