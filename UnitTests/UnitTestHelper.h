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

#ifndef SMOLDOCK_UNITTESTHELPER_H
#define SMOLDOCK_UNITTESTHELPER_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
/*
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/ReactionRunner.h>
#include <GraphMol/ChemReactions/PreprocessRxn.h>
#include <GraphMol/ChemReactions/SanitizeRxn.h>
 */
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>


#include "../Structures/Molecule.h"

namespace SmolDock {

    class UnitTestHelper {
    public:
        UnitTestHelper() = default;

        void populateMol_COCOH(Molecule *mol);

    private:

    };

}

#endif //SMOLDOCK_UNITTESTHELPER_H
