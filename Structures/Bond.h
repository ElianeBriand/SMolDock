//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_BOND_H
#define SMOLDOCK_BOND_H

#include <memory>
#include "Atom.h"

namespace SmolDock {

    class Bond {

    public:
        Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b);


    private:
        std::shared_ptr<Atom> bond_end_a;
        std::shared_ptr<Atom> bond_end_b;
    };

}

#endif //SMOLDOCK_BOND_H
