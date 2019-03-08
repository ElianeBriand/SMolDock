//
// Created by eliane on 31/12/18.
//

#ifndef SMOLDOCK_VINACOMPATIBILITYPOSTPROCESSOR_H
#define SMOLDOCK_VINACOMPATIBILITYPOSTPROCESSOR_H

#include "InputModifierInterface.h"

namespace SmolDock::InputModifier {


    class VinaCompatibility : public InputModifier {

    public:
        VinaCompatibility() = default;


        std::vector<std::tuple<int,int>> deselectRotatableBonds(std::shared_ptr<RDKit::RWMol> rwmol) final;

        void postProcessAtomFromLigand(SmolDock::Atom &atom) final;

        void postProcessAtomFromProtein(SmolDock::Atom &atom, SmolDock::AminoAcid &residue) final;

        ~VinaCompatibility() final = default;

    private:

    };

}


#endif //SMOLDOCK_VINACOMPATIBILITYPOSTPROCESSOR_H
