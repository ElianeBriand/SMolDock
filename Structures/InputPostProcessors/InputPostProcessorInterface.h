//
// Created by eliane on 31/12/18.
//

#ifndef SMOLDOCK_INPUTPOSTPROCESSORINTERFACE_H
#define SMOLDOCK_INPUTPOSTPROCESSORINTERFACE_H

#include <Structures/Atom.h>

namespace SmolDock::InputPostProcessor {

    //! Interface definition for InputPostProcessor
    /*!
     * The purpose of input post processor is mainly to modify ligand and protein atom & residue to simulate the behaviour
     * of a particular software, and get numerically identical results. Though other uses are possible and imaginable.
     *
     * For example, we consider cysteine sulfur atom as a H-bond acceptor, and sometime
     * donor (if not part of a disulfide bridge). However, Vina does not. (Because X-Score does not either). So we correct that
     * in a post processing step.
     *
     * We find that such differences are often rather minimal and do not warrant custom importers.
     *
     */
    class InputPostProcessor {

    public:

        //! Modify the given atom of the ligand
        /*!
         * This fonction will be called for every atom of the ligand
         *
         * It is garantueed that bonds will have been set and propagated to atoms
         * when this function is called. It can use bonded atom for its determinations.
         *
         * \param atom An atom of the ligand
         */
        virtual void processAtomFromLigand(SmolDock::Atom& atom) = 0;

        //! Modify the given atom of the protein
        /*!
         * This fonction will be called for every atom of the protein, with the residue it belongs to
         *
         *
         * \param atom An atom of the protein
         * \param residue The amino acid containing such atom
         */
        virtual void processAtomFromProtein(SmolDock::Atom& atom, SmolDock::AminoAcid& residue) = 0;
    };

}

#endif //SMOLDOCK_INPUTPOSTPROCESSORINTERFACE_H
