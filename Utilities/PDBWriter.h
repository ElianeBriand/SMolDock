//
// Created by eliane on 22/12/18.
//

#ifndef SMOLDOCK_PDBWRITER_H
#define SMOLDOCK_PDBWRITER_H

#include <string>
#include <vector>
#include <Structures/Molecule.h>

namespace SmolDock {

    /*!
     * \brief A PDB file writer for ligands.
     *
     * Not really standard conformant but ChemInfo software will open the generated file without problem.
     *
     */
    class PDBWriter {

    public:
        //! Add a ligand to be exported as PDB
        /*!
         * \param ligand Ligand
        */
        void addLigand(const Molecule& ligand);

        //! Write the added ligand to filename
        /*!
         * \param filename filename to be written. Content is overwritten if existing
         * \sa addLigand()
        */
        bool writePDB(std::string filename);

    private:
        std::vector<Molecule> ligand_to_export;
    };

}


#endif //SMOLDOCK_PDBWRITER_H
