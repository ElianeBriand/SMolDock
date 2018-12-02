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

        enum BondType {
            singlebond,
            doublebond,
            triplebond,
            defaultbond,
            aromatic
        };

        Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b);

        Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b, unsigned int id);

        std::shared_ptr<Atom> getEndA();

        std::shared_ptr<Atom> getEndB();

        void publicizeToAtom();

        void setBondType(Bond::BondType t);

        BondType getBondType() const;

    private:
        std::shared_ptr<Atom> bond_end_a;
        std::shared_ptr<Atom> bond_end_b;

        static unsigned int nextBondID;
        unsigned int BondID;

        BondType bondtype;

    };

}

#endif //SMOLDOCK_BOND_H
