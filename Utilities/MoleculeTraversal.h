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

#ifndef SMOLDOCK_MOLECULETRAVERSAL_H
#define SMOLDOCK_MOLECULETRAVERSAL_H

#include "../Structures/Molecule.h"
#include "../Structures/Atom.h"


namespace SmolDock {


    class MoleculeTraversal {

    public:
        explicit MoleculeTraversal(const Molecule &mol);

        void printTraversal();

    private:
        const Molecule &molecule;

        void traverse(std::shared_ptr<Atom> a);

        std::vector<std::shared_ptr<Atom> > already_visited_atoms;


    };


}

#endif //SMOLDOCK_MOLECULETRAVERSAL_H
