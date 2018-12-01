//
// Created by eliane on 28/11/18.
//

#include <memory>

#include "ConformerRigidDockingEngine.h"
#include "Internals/iConformer.h"


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
/*
        this->rnd_generator.seed(this->random_seed);
        std::uniform_int_distribution<> dis(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());



        std::vector<int> conformer_ids = RDKit::DGeomHelpers::EmbedMultipleConfs(*rwmol,
            this->conformer_num,
            30,
            dis(rnd_generator),
            false);

        std::vector<RDKit::Conformer> rdkit_conformers;
        std::vector<iConformer> conformers;
        rdkit_conformers.reserve(this->conformer_num);
        conformers.reserve(this->conformer_num);

        unsigned int num_atoms = this->orig_ligand->numberOfAtoms();

        for(int i : conformer_ids){
            iConformer conformer;
            RDKit::Conformer rdkit_conformer = this->rwmol->getConformer(i);

            rdkit_conformers.push_back(rdkit_conformer);

            conformer.atoms_vect = std::make_unique< std::vector<iAtom> >(static_cast<size_t>(num_atoms));

            for (auto atom_it = rwmol->beginAtoms(); atom_it != rwmol->endAtoms(); ++atom_it) {
                const RDGeom::Point3D &position = rdkit_conformer.getAtomPos((*atom_it)->getIdx());
            }

            //             const RDGeom::Point3D &position = starting_conformer.getAtomPos((*atom_it)->getIdx());
            //            current_atom->setAtomPosition(std::make_tuple(position.x, position.y, position.z));
            //            atoms.push_back(current_atom);

        }




*/

            return true;
        }

        void ConformerRigidDockingEngine::runDockingEngine() {

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
        }

    }

}