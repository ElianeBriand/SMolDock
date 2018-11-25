//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_MOLECULE_H
#define SMOLDOCK_MOLECULE_H

#include <vector>


#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
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
#include <GraphMol/DistGeomHelpers/Embedder.h>

#include "Structure.h"
#include "Atom.h"
#include "Bond.h"

namespace SmolDock {


    class Molecule : Structure {
        friend class MoleculeTraversal;

        friend class UnitTestHelper;
    public:
        Molecule();



        std::shared_ptr<RDKit::RWMol> getInternalRWMol();

        bool populateFromSMILES(const std::string &smiles);


    private:
        std::vector<std::shared_ptr<Atom> > atoms;
        std::vector<std::shared_ptr<Bond> > bonds;

        std::shared_ptr<RDKit::RWMol> rwmol;
        std::string smiles;

    };

}

#endif //SMOLDOCK_MOLECULE_H
