//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_MOLECULETRAVERSAL_H
#define SMOLDOCK_MOLECULETRAVERSAL_H

#include "../Structures/Molecule.h"
#include "../Structures/Atom.h"


namespace SmolDock {


    class MoleculeTraversal {

    public:
        explicit MoleculeTraversal(const Molecule &mol);

        void printTraversal();

    private:
        const Molecule &molecule;

        void traverse(std::shared_ptr<Atom> a);

        std::vector<std::shared_ptr<Atom> > already_visited_atoms;


    };


}

#endif //SMOLDOCK_MOLECULETRAVERSAL_H
