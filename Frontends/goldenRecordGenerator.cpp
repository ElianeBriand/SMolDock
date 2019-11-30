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

/*
 * This is the golden record generator for the FPGA accelerator project
 *
 * In a nutshell, this generates score and score subcomponent value for known input
 * values, to compare to the FPGA output. This will be especially useful when optimizing
 * the datatypes (impact of replacing doubles with 50-bit floating point, or single precision,
 * or fixed point, in terms on figures difference on the final score)
 *
 */


#include <iostream>
#include <memory>
#include <chrono>

#include <Engines/ScoringFunctions/VinaLike.h>
#include <Engines/ScoringFunctions/VinaLikeRigid.h>
#include <Engines/ScoringFunctions/VinaLikeCommon.h>
#include <Engines/Internals/InternalsUtilityFunctions.h>




#include <Utilities/Version.h>
#include <Utilities/AdvancedErrorHandling.h>
#include <Frontends/FrontendCommon.h>

#include <boost/log/trivial.hpp>

namespace sd = SmolDock;


int main() {

    setupAdvancedErrorHandling();
    sd::setupLogPrinting();


    auto startup_moment = std::chrono::system_clock::now();
    std::time_t startup_time = std::chrono::system_clock::to_time_t(startup_moment);

    BOOST_LOG_TRIVIAL(info) << "SmolDock golden record generator " << sd::getVersionString();
    BOOST_LOG_TRIVIAL(info) << "";
    BOOST_LOG_TRIVIAL(info) << "Data generated on " << std::ctime(&startup_time);
    BOOST_LOG_TRIVIAL(info) << "";
    BOOST_LOG_TRIVIAL(info) << "";


    sd::iConformer conformer1;
    conformer1.num_rotatable_bond = 0;
    conformer1.centroidNormalizingTransform[0] = 0.0;
    conformer1.centroidNormalizingTransform[1] = 0.0;
    conformer1.centroidNormalizingTransform[2] = 0.0;

    conformer1.x.push_back(1.0);
    conformer1.y.push_back(2.0);
    conformer1.z.push_back(0.5);
    conformer1.atomicRadius.push_back(1.7); 
    conformer1.type.push_back(8); // Carbon
    conformer1.variant.push_back(1 << 1);


    BOOST_LOG_TRIVIAL(info) << "Ligand conformation :";

    for(int i = 0; i  < conformer1.x.size(); ++i) {
        BOOST_LOG_TRIVIAL(info) << "ATOM";
        BOOST_LOG_TRIVIAL(info) << "   pos = " << std::fixed<< std::setprecision(5) << conformer1.x.at(i) << ", " << conformer1.y.at(i)  << ", " << conformer1.z.at(i);
        BOOST_LOG_TRIVIAL(info) << "   type = " << std::fixed<< std::setprecision(5) << (unsigned int) conformer1.type.at(i);
        BOOST_LOG_TRIVIAL(info) << "   variant = " << std::fixed<< std::setprecision(5) << conformer1.variant.at(i);
        BOOST_LOG_TRIVIAL(info) << "   radius = " << std::fixed<< std::setprecision(5) << conformer1.atomicRadius.at(i);
    }

    BOOST_LOG_TRIVIAL(info) << "End of ligand conformation";


    sd::iProtein prot1;
    prot1.center_x = 0.0;
    prot1.center_y = 0.0;
    prot1.center_z = 0.0;

    prot1.radius = 10.0;

    prot1.x.push_back(0.0);
    prot1.y.push_back(0.0);
    prot1.z.push_back(0.0);
    prot1.atomicRadius.push_back(1.7);
    prot1.type.push_back(8); // Carbon
    prot1.variant.push_back(1 << 2); 

    BOOST_LOG_TRIVIAL(info) << "Prot conformation :";

    for(int i = 0; i  < prot1.x.size(); ++i) {
        BOOST_LOG_TRIVIAL(info) << "PROT ATOM";
        BOOST_LOG_TRIVIAL(info) << "   pos = " << std::fixed<< std::setprecision(5) << prot1.x.at(i) << ", " << prot1.y.at(i)  << ", " << prot1.z.at(i);
        BOOST_LOG_TRIVIAL(info) << "   type = " << std::fixed<< std::setprecision(5) << (unsigned int)prot1.type.at(i);
        BOOST_LOG_TRIVIAL(info) << "   variant = " << std::fixed<< std::setprecision(5) << prot1.variant.at(i);
        BOOST_LOG_TRIVIAL(info) << "   radius = " << std::fixed<< std::setprecision(5) << prot1.atomicRadius.at(i);
    }

    BOOST_LOG_TRIVIAL(info) << "End of prot conformation";


    sd::iTransform tr1 = sd::iTransformIdentityInit(conformer1.num_rotatable_bond);

    sd::Score::VinaLikeRigid vlr(conformer1,prot1, tr1);

    arma::mat internalRepr =  vlr.getStartingConditions();

    std::vector<std::tuple<std::string,double>>  componentScore = vlr.EvaluateSubcomponents(internalRepr);

    BOOST_LOG_TRIVIAL(info) << "Scoring contribution without weighting :";

    for(int i = 0; i  < componentScore.size(); ++i) {
        BOOST_LOG_TRIVIAL(info) << std::get<std::string>(componentScore.at(i)) << " = " << std::fixed<< std::setprecision(5) << std::get<double>(componentScore.at(i));
    }

    BOOST_LOG_TRIVIAL(info) << "End of record";

    return 0;


}
