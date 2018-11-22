//
// Created by eliane on 11/11/18.
//

#include "Protein.h"
#include "Atom.h"

namespace SmolDock {

    bool Protein::populateFromPDB(const std::string &filename) {


        // Read all non-water atoms/stuff and water in two separate system
        // TODO : do something with the water or switch to one system
        ESBTL::PDB_line_selector_two_systems sel;

        std::vector<ESBTL::Default_system> systems;


        ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems, sel.max_nb_systems());


        //read the pdb file
        if (ESBTL::read_a_pdb_file(filename, sel, builder,
                                   ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >())) {

            if (systems.empty() || systems[0].has_no_model()) {
                std::cerr << "[!] No atoms found in PDB" << std::endl;
                std::exit(-1);
            }

            /**
             * For reference purpose, the hierarchy of the PDB file reading is :
             * System (defined by ESBTL, like input filter in some ways)
             * Model (PDB-file specification : more or less, protein)
             * Chain (Amino acid chain)
             * Residue (=Amino acid)
             * Atom
             *
             * Heteroatoms (Cl, or crystallographic ions) are mixed in at the residue level.
             * We detect those if residueName is not in {set of defined residueName}
             * TODO : Differentiate between legit biological ions (metalloprotein, ...) and crystallographic heavy atoms
             */

            unsigned int total_nb_model = 0, total_nb_chain = 0, total_nb_atom = 0, total_nb_hetatm = 0, total_nb_residue = 0;

            //Iterate over all models in the first system
            for (ESBTL::Default_system::Models_iterator
                         it_model = systems[0].models_begin();
                 it_model != systems[0].models_end();
                 ++it_model) {
                const ESBTL::Default_system::Model &model = *it_model;
                total_nb_model++;

                // Iterating over chains
                for (ESBTL::Default_system::Model::Chains_const_iterator it_ch = model.chains_begin();
                     it_ch != model.chains_end(); ++it_ch) {
                    total_nb_chain++;
                    // Iterating over residues in the chain
                    for (ESBTL::Default_system::Chain::Residues_const_iterator it_res = it_ch->residues_begin();
                         it_res != it_ch->residues_end(); ++it_res) {
                        total_nb_residue++;
                        if (stringToResType(it_res->residue_name()) == AminoAcid::AAType::unknown) {
                            // We have an heteroatom
                            this->heteroatoms.emplace_back(std::shared_ptr<Atom>(new Atom(it_res->residue_name())));
                            total_nb_hetatm++;
                            total_nb_atom++;
                            continue;
                        }
                        // It's a normal residue
                        auto current_residue = this->aminoacids.emplace_back(
                                std::shared_ptr<AminoAcid>(new AminoAcid(it_res->residue_name()))
                        );
                        // Iterate over atoms in the residue
                        for (ESBTL::Default_system::Residue::Atoms_const_iterator it_atm = it_res->atoms_begin();
                             it_atm != it_res->atoms_end(); ++it_atm) {
                            assert(it_atm->is_hetatm() ==
                                   false); // We expect heteroatoms to have been taken care of previously
                            total_nb_atom++;
                            auto current_atom = current_residue->atoms.emplace_back(
                                    std::shared_ptr<Atom>(new Atom(it_atm->atom_name(), true,
                                                                   stringToResType(it_res->residue_name())))
                            );
                            current_atom->setAtomPosition(std::make_tuple(it_atm->x(), it_atm->y(), it_atm->z()));
                        }
                    }
                }
                //Consider only the first model of the first system
                std::cout << filename << " : loaded " << total_nb_model << " models, "
                          << total_nb_chain << " chains, " << total_nb_residue << " residues, "
                          << total_nb_atom << " atoms (" << total_nb_hetatm << " heteroatoms)." << std::endl;
            }
        } else
            return false;

        return true;
    }


}