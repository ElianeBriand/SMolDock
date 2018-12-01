//
// Created by eliane on 28/11/18.
//

#ifndef SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H
#define SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H

#include <GraphMol/RWMol.h>
#include <random>


#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"


namespace SmolDock::Engine {

/*
 * ConformerRigidDockingEngine
 *
 * The general idea is generating lots of conformer, then docking them with constant bond angle
 * and distance ("rigid" docking). The theory being that this may be faster than real docking
 * and still produce a similar affinity ranking.
 *
 * )
 *
 */
    //! Generates a lot of conformer, then runs rigid docking on them.
    /*!
     * The general idea is generating lots of conformer, then docking them with constant bond angle
     * and distance ("rigid" docking). The theory being that this may be faster than real docking
     * and still produce a similar affinity ranking. (In practice, I do not know if it would work :
     * hence this to try it
     */
    class ConformerRigidDockingEngine : public AbstractDockingEngine {

    public:
        ConformerRigidDockingEngine(unsigned int conformer_num);

        ///// Parameters /////////////
        bool setDockingBox(DockingBoxSetting setting) final;

        bool setProtein(Protein *p) final;

        bool setLigand(Molecule *m) final;

        void setRandomSeed(int seed) final;


        ///// Actions /////////////
        bool setupDockingEngine() final;

        void runDockingEngine() final;

        ///// Results /////////////
        std::shared_ptr<DockingResult> getDockingResult() final;


    private:

        unsigned int conformer_num;

        Protein *orig_protein;
        Molecule *orig_ligand;

        int random_seed = 1;
        std::mt19937 rnd_generator;

        std::shared_ptr<RDKit::RWMol> rwmol;
    };

}


#endif //SMOLDOCK_CONFORMERRIGIDDOCKINGENGINE_H
