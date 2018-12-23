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
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef SMOLDOCK_BASICSCORINGFUNCTION_H
#define SMOLDOCK_BASICSCORINGFUNCTION_H


#include <Engines/Internals/iConformer.h>
#include <Engines/Internals/iProtein.h>
#include <Engines/Internals/iTransform.h>


namespace SmolDock {

    /*! \namespace SmolDock::Score Free-standing docking score functions */
    namespace Score {

        //! \fn Score a protein-ligand configuration. Variant : Basic
        /*!
         * Process both the intra- and inter-molecular part of the docking score.
         * The Basic variant is a  simple scoring function destined for
         * general-purpose affinity ranking.
         *
         * \param conformer Ligand conformation & position to evaluate
         * \param transform Transformation to apply to the ligand
         * \param protein Protein conformation to evaluate
         * \return The docking score
         * \sa
        */
        double basic_scoring_func(iConformer& conformer, iTransform& transform, iProtein& protein);


    }

}


#endif //SMOLDOCK_BASICSCORINGFUNCTION_H