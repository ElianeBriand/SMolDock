//
// Created by eliane on 16/03/19.
//

#ifndef SMOLDOCK_SMARTSMATCHER_H
#define SMOLDOCK_SMARTSMATCHER_H

#include <string>
#include <memory>

#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>

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
