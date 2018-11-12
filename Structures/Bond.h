//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_BOND_H
#define SMOLDOCK_BOND_H

#include <memory>
#include "Atom.h"

namespace SmolDock {

    class Bond : public std::enable_shared_from_this<Bond> {

    public:
        Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b);

        std::shared_ptr<Atom> getEndA();

        std::shared_ptr<Atom> getEndB();

        void publicizeToAtom();
    private:
        std::shared_ptr<Atom> bond_end_a;
        std::shared_ptr<Atom> bond_end_b;
    };

}

#endif //SMOLDOCK_BOND_H
