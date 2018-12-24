//
// Created by eliane on 24/12/18.
//

#ifndef SMOLDOCK_INTERMEDIATECONFORMERCOLLECTOR_H
#define SMOLDOCK_INTERMEDIATECONFORMERCOLLECTOR_H

#include <Structures/Molecule.h>
#include "PDBWriter.h"

namespace SmolDock {

class IntermediateConformerCollector {
public:
    IntermediateConformerCollector(Molecule* mol, PDBWriter* pdbwriter);
    void addiConformer(iConformer conformer);

private:
    Molecule* molecule = nullptr;
    PDBWriter* writer = nullptr;


};

}


#endif //SMOLDOCK_INTERMEDIATECONFORMERCOLLECTOR_H
