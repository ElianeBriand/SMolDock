//
// Created by eliane on 31/12/18.
//

#ifndef SMOLDOCK_VINACOMPATIBILITYPOSTPROCESSOR_H
#define SMOLDOCK_VINACOMPATIBILITYPOSTPROCESSOR_H

#include "InputPostProcessorInterface.h"

namespace SmolDock::InputPostProcessor {


    class VinaCompatibilityPostProcessor : public InputPostProcessor{

    public:
        VinaCompatibilityPostProcessor() = default;

        void processAtomFromLigand(SmolDock::Atom& atom) final;
        void processAtomFromProtein(SmolDock::Atom& atom, SmolDock::AminoAcid& residue) final;
    private:

    };

}


#endif //SMOLDOCK_VINACOMPATIBILITYPOSTPROCESSOR_H
