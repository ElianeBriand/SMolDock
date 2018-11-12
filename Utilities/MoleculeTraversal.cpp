//
// Created by eliane on 11/11/18.
//

#include <iostream>
#include "MoleculeTraversal.h"

using namespace std;

namespace SmolDock {


    MoleculeTraversal::MoleculeTraversal(const Molecule &mol) : molecule(mol) {

    }

    void MoleculeTraversal::traverse(shared_ptr<Atom> a) {
        cout << a->getAtomID() << " : " << a->getTypeString() << " (" << a->bonds.size() << " bonds)" << endl;
        for (auto bond: a->bonds) {
            auto b = bond.lock();
            if (b == nullptr) {
                cout << "MoleculeTraversal::traverse : Cant lock weak ptr (?)" << endl;
            }
            auto endA = b->getEndA();
            auto endB = b->getEndB();
            auto nonself_end = (endA == a) ? endB : endA; // Which end of the bond is not this atom

            for (auto already_visited_atom: already_visited_atoms) {
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