//
// Created by eliane on 23/12/18.
//

#ifndef SMOLDOCK_PDBLIGANDUTILS_H
#define SMOLDOCK_PDBLIGANDUTILS_H

#include <memory>

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

namespace SmolDock {

    void AssignBondOrderFromTemplateSMILES(std::shared_ptr<RDKit::RWMol> mol, const std::string &smiles);

}
#endif //SMOLDOCK_PDBLIGANDUTILS_H
