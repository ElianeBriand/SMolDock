//
// Created by eliane on 14/11/18.
//

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

        void populateMol_COCOH(std::shared_ptr<Molecule> mol);

        std::tuple<unsigned long, unsigned long> getNumberOfAtomAndBonds(std::shared_ptr<Molecule> mol);

        RDKit::RWMol *getRDKitRWMol(std::shared_ptr<Molecule> mol);

    private:

    };

}

#endif //SMOLDOCK_UNITTESTHELPER_H
