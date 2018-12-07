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
