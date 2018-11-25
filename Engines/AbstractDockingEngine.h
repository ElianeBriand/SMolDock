//
// Created by eliane on 25/11/18.
//

#ifndef SMOLDOCK_ABSTRACTDOCKINGENGINE_H
#define SMOLDOCK_ABSTRACTDOCKINGENGINE_H

#include "../Structures/Protein.h"
#include "../Structures/Results/DockingResult.h"
#include "../Structures/Molecule.h"

namespace SmolDock {

    class AbstractDockingEngine {
    public:

        enum class DockingBoxSetting {
            everything,
            solventExposed
        };

        enum class DockingBoxShape {
            sphere,
            cube
        };

        virtual bool setProtein(const Protein *p) = 0;

        virtual bool setMolecule(const Molecule *m) = 0;

        virtual bool setDockingBox(DockingBoxSetting setting) = 0;
        // virtual bool setDockingBox(std::tuple<double,double,double> center, double r, DockingBoxShape shape = DockingBoxShape::sphere) = 0;



        virtual bool setupDockingEngine() = 0;

        virtual void runDockingEngine() = 0;

        virtual std::shared_ptr<DockingResult> getDockingResult() = 0;
    };

}

#endif //SMOLDOCK_ABSTRACTDOCKINGENGINE_H
