//
// Created by eliane on 11/11/18.
//

#include "Bond.h"

namespace SmolDock {


    Bond::Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b) {
        this->bond_end_a = atom_a;
        this->bond_end_b = atom_b;
    }
}
