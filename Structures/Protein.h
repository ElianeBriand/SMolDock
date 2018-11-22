//
// Created by eliane on 11/11/18.
//

#ifndef SMOLDOCK_PROTEIN_H
#define SMOLDOCK_PROTEIN_H


#include <vector>
#include <memory>

#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>

#include "Structure.h"
#include "AminoAcid.h"

namespace SmolDock {


    class Protein : Structure {

    public:
        Protein() = default;

        bool populateFromPDB(const std::string &filename);

    private:
        std::vector<std::shared_ptr<AminoAcid> > aminoacids;
        std::vector<std::shared_ptr<Atom> > heteroatoms;

    };

}

#endif //SMOLDOCK_PROTEIN_H
