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

#include <iostream>
#include "MoleculeTraversal.h"

using namespace std;

namespace SmolDock {


    MoleculeTraversal::MoleculeTraversal(const Molecule &mol) : molecule(mol) {

    }

    void MoleculeTraversal::traverse(shared_ptr<Atom> a) {
        cout << a->getAtomID() << " : " << a->getTypeAsString() << " (" << a->bonds.size() << " bonds)" << endl;
        for (const auto &bond: a->bonds) {
            auto b = bond.lock();
            if (b == nullptr) {
                cout << "MoleculeTraversal::traverse : Cant lock weak ptr (?)" << endl;
            }
            auto endA = b->getEndA();
            auto endB = b->getEndB();
            auto nonself_end = (endA == a) ? endB : endA; // Which end of the bond is not this atom

            for (const auto &already_visited_atom: already_visited_atoms) {
                if (nonself_end == already_visited_atom) // already visited
                {
                    goto nextbond; // end traversal, go to next iteration
                }
            }
            already_visited_atoms.push_back(a);
            traverse(nonself_end); // Continue deeper
            nextbond:;
        }
    }

    void MoleculeTraversal::printTraversal() {
        cout << "=== MoleculeTraversal ===" << endl;
        size_t n_of_atom = molecule.atoms.size();

        if (n_of_atom == 0) {
            cout << "Empty molecule" << endl;
            cout << "=== END ===" << endl;
            return;

        }

        auto first_atom = molecule.atoms[0];
        traverse(first_atom);

        cout << "=== END ===" << endl;
    }
}