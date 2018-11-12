//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_MOLECULE_H
#define SMOLDOCK_MOLECULE_H

#include <vector>
#include "Structure.h"
#include "Atom.h"
#include "Bond.h"

namespace SmolDock {

    class MoleculeTraversal;

    class Molecule : Structure {
        friend MoleculeTraversal;
    public:
        Molecule();


        /*** NOT FOR USE IN ACTUAL CODE ****/
        void _dev_populateSampleMolecule();


    private:
        std::vector<std::shared_ptr<Atom> > atoms;
        std::vector<std::shared_ptr<Bond> > bonds;
    };

}

#endif //SMOLDOCK_MOLECULE_H
