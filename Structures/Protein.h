//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_PROTEIN_H
#define SMOLDOCK_PROTEIN_H


#include <vector>
#include <memory>

#include "Structure.h"
#include "AminoAcid.h"

namespace SmolDock {


    class Protein : Structure {

    private:
        std::vector<std::shared_ptr<AminoAcid> > aminoacids;

    };

}

#endif //SMOLDOCK_PROTEIN_H
