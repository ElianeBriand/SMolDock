//
// Created by eliane on 25/11/18.
//

#ifndef SMOLDOCK_MDSTYLEDOCKINGENGINE_H
#define SMOLDOCK_MDSTYLEDOCKINGENGINE_H

#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"

namespace SmolDock::Engine {

    class MDStyleDockingEngine : AbstractDockingEngine {


    public:
        MDStyleDockingEngine() = default;


        bool setProtein(Protein *p) final;

        bool setLigand(Molecule *m) final;


        bool setDockingBox(DockingBoxSetting setting) final;

        bool setupDockingEngine() final;

        void runDockingEngine() final;

        std::shared_ptr<DockingResult> getDockingResult() final;


    private:

        Protein *protein = nullptr;
        Molecule *ligand = nullptr;

    };

}

#endif //SMOLDOCK_MDSTYLEDOCKINGENGINE_H
