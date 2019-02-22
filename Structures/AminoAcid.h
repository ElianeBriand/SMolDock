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

#ifndef SMOLDOCK_AMINOACID_H
#define SMOLDOCK_AMINOACID_H


#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <memory>

#include "Engines/Internals/iProtein.h"


namespace SmolDock {


    class Atom;

    enum class PDBResidueVariantAssignationType;

    class AminoAcid {

        friend class Protein;

        friend void
        assignVariantFlagsForResidueAtom(AminoAcid &residue, PDBResidueVariantAssignationType assignation_type);

    public:

        enum class AAType {
            alanine,
            arginine,
            asparagine,
            aspartate,
            cysteine,
            glutamate,
            glutamine,
            glycine,
            histidine,
            isoleucine,
            leucine,
            lysine,
            methionine,
            phenylalanine,
            proline,
            serine,
            threonine,
            tryptophan,
            tyrosine,
            valine,
            selenocysteine,
            pyrrolysine,
            heteroatom, /* Used for convenience */
            unknown
        };

        /// Internal type <=> string conversion //////
        friend std::string resTypeToString(AminoAcid::AAType t);

        friend AminoAcid::AAType stringToResType(const std::string &shorthand_or_name);


        explicit AminoAcid(const std::string &AA3LettersShorthand);

        explicit AminoAcid(const AminoAcid::AAType &type);


        unsigned int getAAId() const;

        void setAAId(unsigned int id);

        bool filliProtein(iProtein &prot, bool skipHydrogen = true);

        AAType getType() const;

        void setType(AAType t);

    protected:
        static std::set<std::tuple<AminoAcid::AAType, std::string, std::string> > AAShorthandSet;

        std::array<double, 3> centroid;
        double maxDistanceFromCentroid;


    private:
        AAType type;

        unsigned int AAId;

        std::string shorthand;
        std::string fullName;

        std::vector<std::shared_ptr<Atom> > atoms;


    };

    std::string resTypeToString(AminoAcid::AAType t);

    AminoAcid::AAType stringToResType(const std::string &shorthand_or_name);

}


#endif //SMOLDOCK_AMINOACID_H
