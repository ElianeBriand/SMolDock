//
// Created by eliane on 16/03/19.
//

#ifndef SMOLDOCK_SMARTSMATCHER_H
#define SMOLDOCK_SMARTSMATCHER_H

#include <string>
#include <memory>

#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/FileParsers/MolSupplier.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

namespace SmolDock {

    class SMARTSMatcher {
    public:
        SMARTSMatcher(const std::string& smarts_);

        bool matchesSMILES(const std::string& smiles_);

    private:
        std::string smarts;

        std::shared_ptr<RDKit::RWMol> pattern;


    };
}


#endif //SMOLDOCK_SMARTSMATCHER_H
