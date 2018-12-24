//
// Created by eliane on 24/12/18.
//

#include "IntermediateConformerCollector.h"


namespace SmolDock {




    void IntermediateConformerCollector::addiConformer(iConformer conformer) {
        Molecule intermediate_res = this->molecule->deepcopy();
        intermediate_res.updateAtomPositionsFromiConformer(conformer);
        writer->addLigand(intermediate_res);
    }

    IntermediateConformerCollector::IntermediateConformerCollector(Molecule *mol, PDBWriter *pdbwriter):
    molecule(mol),writer(pdbwriter){

    }
}