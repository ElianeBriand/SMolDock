//
// Created by eliane on 25/11/18.
//

#ifndef SMOLDOCK_MDSTYLEDOCKINGENGINE_H
#define SMOLDOCK_MDSTYLEDOCKINGENGINE_H

#include "AbstractDockingEngine.h"
#include "../Structures/Protein.h"
#include "../Structures/Molecule.h"

namespace SmolDock {

    class MDStyleDockingEngine : AbstractDockingEngine {


    public:
        MDStyleDockingEngine() = default;


        bool setProtein(const Protein *p) final;

        bool setMolecule(const Molecule *m) final;


        bool setDockingBox(DockingBoxSetting setting) final;

        bool setupDockingEngine() final;

        void runDockingEngine() final;

        std::shared_ptr<DockingResult> getDockingResult() final;


    private:

        const Protein *protein;
        const Molecule *molecule;

    };

}

#endif //SMOLDOCK_MDSTYLEDOCKINGENGINE_H
