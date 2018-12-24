/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "Bond.h"

namespace SmolDock {

    unsigned int Bond::nextBondID = 0;


    Bond::Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b) {
        this->bond_end_a = std::move(atom_a);
        this->bond_end_b = std::move(atom_b);
    }

    Bond::Bond(std::shared_ptr<Atom> atom_a, std::shared_ptr<Atom> atom_b, unsigned int id) {
        this->bond_end_a = std::move(atom_a);
        this->bond_end_b = std::move(atom_b);

        BondID = id;
        bondtype = defaultbond;
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

    void Bond::setBondType(const Bond::BondType t) {
        bondtype = t;
    }

    Bond::BondType Bond::getBondType() const {
        return bondtype;
    }

    unsigned int Bond::getBondID() {
        return this->BondID;
    }


}