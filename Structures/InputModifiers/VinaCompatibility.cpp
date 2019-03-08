//
// Created by eliane on 31/12/18.
//

#include "VinaCompatibility.h"


namespace SmolDock::InputModifier {


    void VinaCompatibility::postProcessAtomFromLigand(SmolDock::Atom &atom) {
        // We have no special ligand-related work-arounds

    }

    void VinaCompatibility::postProcessAtomFromProtein(SmolDock::Atom &atom, SmolDock::AminoAcid &residue) {

        // Vina follow X-score convention and do not consider sulfur as acceptor or donor
        // (We add this at the PDB ATOM name variant assignation stage : SmolDock::assignPropertiesForResidueAtom)
        if ((residue.getType() == AminoAcid::AAType::cysteine)
            && (atom.getAtomType() == Atom::AtomType::sulfur)) {
            auto atomvar = static_cast<unsigned int>(atom.getAtomVariant());
            unsigned int newatomvar = atomvar & ~(static_cast<unsigned int>(Atom::AtomVariant::hydrogenDonor));
            newatomvar &= ~(static_cast<unsigned int>(Atom::AtomVariant::hydrogenAcceptor));
            atom.setAtomVariant(static_cast<Atom::AtomVariant>(newatomvar));
        }
    }

    std::vector<std::tuple<int,int>> VinaCompatibility::deselectRotatableBonds(std::shared_ptr<RDKit::RWMol> rwmol) {
        std::vector<std::tuple<int,int>> empty;
        return empty;
    }
}