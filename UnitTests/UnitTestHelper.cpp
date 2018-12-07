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

#include "UnitTestHelper.h"

namespace SmolDock {

    /* Molecule created : C-O-C-OH */
    void UnitTestHelper::populateMol_COCOH(Molecule *mol) {
        auto atom1 = std::make_shared<Atom>(Atom::AtomType::carbon);
        mol->atoms.push_back(atom1);

        auto atom1_h1 = std::make_shared<Atom>(Atom::AtomType::hydrogen);
        auto bond_atom1_h1 = std::make_shared<Bond>(atom1, atom1_h1);
        mol->atoms.push_back(atom1_h1);
        mol->bonds.push_back(bond_atom1_h1);


        auto atom1_h2 = std::make_shared<Atom>(Atom::AtomType::hydrogen);
        auto bond_atom1_h2 = std::make_shared<Bond>(atom1, atom1_h2);
        mol->atoms.push_back(atom1_h2);
        mol->bonds.push_back(bond_atom1_h2);

        auto atom1_h3 = std::make_shared<Atom>(Atom::AtomType::hydrogen);
        auto bond_atom1_h3 = std::make_shared<Bond>(atom1, atom1_h3);
        mol->atoms.push_back(atom1_h3);
        mol->bonds.push_back(bond_atom1_h3);

        auto atom2 = std::make_shared<Atom>(Atom::AtomType::oxygen);
        auto bond_atom1_atom2 = std::make_shared<Bond>(atom1, atom2);
        mol->atoms.push_back(atom2);
        mol->bonds.push_back(bond_atom1_atom2);

        auto atom3 = std::make_shared<Atom>(Atom::AtomType::carbon);
        auto bond_atom2_atom3 = std::make_shared<Bond>(atom2, atom3);
        mol->atoms.push_back(atom3);
        mol->bonds.push_back(bond_atom2_atom3);

        auto atom3_h1 = std::make_shared<Atom>(Atom::AtomType::hydrogen);
        auto bond_atom3_h1 = std::make_shared<Bond>(atom3, atom3_h1);
        mol->atoms.push_back(atom3_h1);
        mol->bonds.push_back(bond_atom3_h1);

        auto atom3_h2 = std::make_shared<Atom>(Atom::AtomType::hydrogen);
        auto bond_atom3_h2 = std::make_shared<Bond>(atom3, atom3_h2);
        mol->atoms.push_back(atom3_h2);
        mol->bonds.push_back(bond_atom3_h2);

        auto atom4 = std::make_shared<Atom>(Atom::AtomType::oxygen);
        auto bond_atom3_atom4 = std::make_shared<Bond>(atom3, atom4);
        mol->atoms.push_back(atom4);
        mol->bonds.push_back(bond_atom3_atom4);

        auto atom5 = std::make_shared<Atom>(Atom::AtomType::hydrogen);
        auto bond_atom4_atom5 = std::make_shared<Bond>(atom4, atom5);
        mol->atoms.push_back(atom5);
        mol->bonds.push_back(bond_atom4_atom5);


        for (const auto &bond: mol->bonds) {
            bond->publicizeToAtom();
        }
    }


}