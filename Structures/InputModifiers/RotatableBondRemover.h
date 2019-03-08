//
// Created by briand on 3/8/19.
//

#ifndef SMOLDOCK_ROTATABLEBONDREMOVER_H
#define SMOLDOCK_ROTATABLEBONDREMOVER_H

#include "InputModifierInterface.h"


namespace SmolDock::InputModifier {


    class RotatableBondRemover : public InputModifier {

    public:
        RotatableBondRemover(std::string SMARTS_pattern_);


        std::vector<std::tuple<int,int>> deselectRotatableBonds(std::shared_ptr<RDKit::RWMol> rwmol) final;

        void postProcessAtomFromLigand(SmolDock::Atom &atom) final;

        void postProcessAtomFromProtein(SmolDock::Atom &atom, SmolDock::AminoAcid &residue) final;

        ~RotatableBondRemover() final = default;

    private:
        std::string SMARTS_pattern;
        std::shared_ptr<RDKit::RWMol> flaggingPattern;
        int idxOfFirstAtomMappedAtom = -1;
        int idxOfSecondAtomMappedAtom = -1;
        bool invalid = false;

    };

}



#endif //SMOLDOCK_ROTATABLEBONDREMOVER_H
