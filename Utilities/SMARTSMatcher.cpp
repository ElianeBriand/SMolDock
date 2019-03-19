//
// Created by eliane on 16/03/19.
//

#include "SMARTSMatcher.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

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


