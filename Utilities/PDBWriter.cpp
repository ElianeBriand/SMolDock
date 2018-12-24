//
// Created by eliane on 22/12/18.
//

#include "PDBWriter.h"

#include <iomanip>

#include <boost/log/trivial.hpp>


#include <ESBTL/default.h>
#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/compressed_ifstream.h>


namespace SmolDock {


    bool PDBWriter::writePDB(std::string filename) {

        ESBTL::Default_system system(0, "defsyst");
        std::ostringstream filecontent;

        filecontent << "AUTHOR    SmolDock" << std::endl;

        int current_model = 1;
        for(const auto& ligand : this->ligand_to_export)
        {
                ESBTL::Molecular_model<ESBTL::Default_system> model(0, system);
                ESBTL::Molecular_chain<ESBTL::Default_system> chain('A', model);
                ESBTL::Molecular_residue<ESBTL::Default_system> residue(ligand.getResidueName(), 0, 'A', chain);

                filecontent << "MODEL      " << std::setw(3) << current_model <<std::endl;
                current_model++;

                int current_atom_serial = 0;
                for(const auto& atom_shptr: ligand.atoms)
                {
                    auto pos = atom_shptr->getAtomPosition();
                    ESBTL::Molecular_atom<ESBTL::Default_system, ESBTL::Point_3> atom1(std::get<0>(pos),std::get<1>(pos), std::get<2>(pos));
                    atom1.is_hetatm() = true;
                    atom1.atom_serial_number() = current_atom_serial;
                    current_atom_serial++;
                    atom1.atom_name() = atomTypeToSymbolString(atom_shptr->getAtomType());
                    atom1.alternate_location() = ' ';
                    atom1.occupancy() = 1.0;
                    atom1.temperature_factor() = 0.0;
                    atom1.element() =  atomTypeToSymbolString(atom_shptr->getAtomType());
                    atom1.charge() = atom_shptr->getCharge();
                    atom1.residue_ = &residue;
                    filecontent << ESBTL::PDB::get_atom_pdb_format(atom1) << std::endl;
                }
                filecontent << "TER   " << std::setw(5) << current_atom_serial << "      LIG A   0A" << std::endl;

                filecontent << "ENDMDL" << std::endl;

        }

        /*
        std::cout << "##### PDB #####" << std::endl;
       std::cout << filecontent.str() << std::endl;
        std::cout << "##### PDB #####" << std::endl;
        */

        std::ofstream pdbfile;
        pdbfile.open(filename);
        pdbfile << filecontent.str() << std::endl;
        pdbfile.close();

        BOOST_LOG_TRIVIAL(info) << "Wrote PDB file :  " << filename;

        return true;
    }

    void PDBWriter::addLigand(const Molecule& ligand) {
            this->ligand_to_export.push_back(ligand);
    }

}
