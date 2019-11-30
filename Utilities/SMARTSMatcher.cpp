//
// Created by eliane on 16/03/19.
//

#include "SMARTSMatcher.h"

#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

namespace SmolDock {
    SMARTSMatcher::SMARTSMatcher(const std::string &smarts_):
    smarts(smarts_),
    pattern(RDKit::SmartsToMol(smarts_)){
    }

    bool SMARTSMatcher::matchesSMILES(const std::string &smiles_) {
        std::shared_ptr<RDKit::ROMol> mol1(RDKit::SmilesToMol(smiles_));

        RDKit::MatchVectType res;
        return RDKit::SubstructMatch(*mol1 , *this->pattern , res );

    }
}


