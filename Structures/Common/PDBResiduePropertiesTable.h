//
// Created by eliane on 29/12/18.
//

#ifndef SMOLDOCK_PDBRESIDUEVARIANTTABLE_H
#define SMOLDOCK_PDBRESIDUEVARIANTTABLE_H

#include <set>
#include <tuple>
#include <vector>
#include <string>

#include <Structures/AminoAcid.h>
#include <Structures/Atom.h>

namespace SmolDock {

    enum class PDBResidueVariantAssignationType {
        GeneralPurpose

    };

    // See the cpp file
    extern std::set<std::tuple<AminoAcid::AAType, std::vector<std::tuple<std::string, Atom::AtomVariant, double> > > > ResidueAtomPropertiesLookupTable_General;


}


#endif //SMOLDOCK_PDBRESIDUEVARIANTTABLE_H
