/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SmolDock is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SmolDock.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <memory>
#include <boost/log/trivial.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include <Utilities/TimingsLog.h>

#include "Protein.h"
#include "Atom.h"

#include "Structures/Common/VariantFlagAssignation.h"
#include "Structures/Common/PDBResidueVariantTable.h"

namespace SmolDock {

    bool Protein::populateFromPDB(const std::string &filename, std::vector<InputPostProcessor::InputPostProcessor*> postProcessors) {



        // Read all non-water atoms/stuff and water in two separate system
        // TODO : do something with the water or switch to one system
        ESBTL::PDB_line_selector_two_systems sel;

        std::vector<ESBTL::Default_system> systems;


        ESBTL::All_atom_system_builder<ESBTL::Default_system> builder(systems, sel.max_nb_systems());


        if (ESBTL::read_a_pdb_file(filename, sel, builder,
                                   ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >())) {

            if (systems.empty() || systems[0].has_no_model()) {
                BOOST_LOG_TRIVIAL(error) << "No atoms found in protein PDB";
                std::exit(-1);
            }

            /*
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

            using namespace boost::accumulators;
            accumulator_set<double, stats<tag::mean, tag::count, tag::min, tag::max > > acc_x, acc_y, acc_z;


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
                            this->heteroatoms.emplace_back(std::make_shared<Atom>(it_res->residue_name()));
                            total_nb_hetatm++;
                            total_nb_atom++;
                            continue;
                        }

                        // else : It's a normal residue
                        auto current_residue = this->aminoacids.emplace_back(
                                std::shared_ptr<AminoAcid>(new AminoAcid(it_res->residue_name()))
                        );

                        current_residue->setAAId(total_nb_residue);

                        // Iterate over atoms in the residue
                        for (ESBTL::Default_system::Residue::Atoms_const_iterator it_atm = it_res->atoms_begin();
                             it_atm != it_res->atoms_end(); ++it_atm) {

                            assert(it_atm->is_hetatm() == 0); // We expect heteroatoms to have been taken care of previously


                            total_nb_atom++;

                            auto current_atom = current_residue->atoms.emplace_back(
                                    std::shared_ptr<Atom>(new Atom(it_atm->atom_name(), true,
                                                                   stringToResType(it_res->residue_name())))
                            );

                            current_atom->setAtomPosition(std::make_tuple(it_atm->x(), it_atm->y(), it_atm->z()));


                            acc_x(it_atm->x());
                            acc_y(it_atm->y());
                            acc_z(it_atm->z());



                        }

                        assignVariantFlagsForResidueAtom(*current_residue, PDBResidueVariantAssignationType::GeneralPurpose);

                    }
                }


                // Post processing
                for(InputPostProcessor::InputPostProcessor* pprocessor : postProcessors )
                {
                    for(auto& residue: aminoacids)
                    {
                        for(auto& atom: residue->atoms)
                        {
                            pprocessor->processAtomFromProtein(*atom,*residue);
                        }
                    }
                }


                this->center_x = mean(acc_x);
                this->center_y = mean(acc_y);
                this->center_z = mean(acc_z);

                this->min_x = min(acc_x);
                this->min_y = min(acc_y);
                this->min_z = min(acc_z);

                this->max_x = max(acc_x);
                this->max_y = max(acc_y);
                this->max_z = max(acc_z);


                BOOST_LOG_TRIVIAL(info) << "Loaded protein : " << filename << "\n    --> " << total_nb_model << " models, "
                          << total_nb_chain << " chains, " << total_nb_residue << " residues, "
                          << total_nb_atom << " atoms (" << total_nb_hetatm << " heteroatoms)";

#ifdef SMOLDOCK_VERBOSE_DEBUG
                BOOST_LOG_TRIVIAL(debug) << "Protein centered on : ("
                          << this->center_x << ", " << this->center_y << ", " <<   this->center_z << ") [n=" <<
                          count(acc_x) << "]";
#endif

            }
        } else {
            BOOST_LOG_TRIVIAL(error) << "Could not load protein PDB file: " << filename ;
            return false;
        }


        // Post processing steps to add variant flags
/*
        for(auto& residue : this->aminoacids)
        {
            //residue->atoms
        }
*/


        return true;
    }

    iProtein Protein::getiProtein() const {
        /*
         *     struct iProtein {
        //! Pseudo-center of protein as a mean of each atom coordinate
        double center_x,center_y,center_z;
        unsigned long size;
        std::vector<double> x,y,z;
        std::vector<unsigned char> type;
        //! Map the AAId of residues to the corresponding position of the atoms in the vector
        std::map<unsigned int, std::tuple<unsigned long,unsigned long> > AAId_to_AtomPositionInVect;
    };
         */
        iProtein prot;

        prot.center_x = this->center_x;
        prot.center_y = this->center_y;
        prot.center_z = this->center_z;


        for(auto& residue: this->aminoacids)
        {
            unsigned long size_before = prot.x.size();
            residue->filliProtein(prot, true /* Skip hydrogen */);
            unsigned long size_after = prot.x.size();

            // Visual demo :
            //
            // index
            //
            // size 0
            //
            // index    0   1
            //         [A] [A]   For residue A : <begin_idx, end_idx> = <0,1> = <size_before, size_after -1>
            // size 0 > 1 > 2
            //
            // index    0   1   2   3   4
            //         [A] [A] [B] [B] [B] For residue B : <begin_idx, end_idx> = <2,4> = <size_before, size_after -1>
            // size 0 > 1 > 2 > 3 > 4 > 5
            //
            // The size_before gives us the index of the first atom that belongs to this residue
            // The size_after lags behind by one.
            // This works even if there is only one atom added (then begin_idx = end_idx = size_before = size_after-1)
            // However there is an edge case : when the residue has no atom, this is incorrect, so :
            assert(size_after > size_before); // At least one atom added

            prot.AAId_to_AtomPositionInVect[residue->getAAId()] = std::make_tuple(size_before, size_after - 1);
        }

        return prot;
    }


}