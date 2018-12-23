//
// Created by eliane on 23/12/18.
//

#ifndef SMOLDOCK_PDBLIGANDUTILS_H
#define SMOLDOCK_PDBLIGANDUTILS_H

#include <memory>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace SmolDock {

void AssignBondOrderFromTemplateSMILES(std::shared_ptr<RDKit::RWMol> mol,const std::string& smiles);

}
#endif //SMOLDOCK_PDBLIGANDUTILS_H
