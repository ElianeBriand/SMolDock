//
// Created by eliane on 11/11/18.
//

#include "Bond.h"

namespace SmolDock {


    Bond::Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b) {
        this->bond_end_a = atom_a;
        this->bond_end_b = atom_b;


    }

    std::shared_ptr<Atom> Bond::getEndA() {
        return bond_end_a;
    }

    std::shared_ptr<Atom> Bond::getEndB() {
        return bond_end_b;
    }

    void Bond::publicizeToAtom() {
        // Publicize the bond to the atoms involved
        // (allow for easy atom based molecule traversal)
        bond_end_a->bonds.push_back(shared_from_this());
        bond_end_b->bonds.push_back(shared_from_this());
    }

}