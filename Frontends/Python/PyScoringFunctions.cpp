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

#include <boost/python.hpp>
#include <boost/numpy.hpp>

#include "PyScoringFunctions.h"

#include <Engines/ScoringFunctions/VinaLikeScoringFunction.h>


namespace p = boost::python;
namespace np = boost::numpy;

namespace sd = SmolDock;

namespace SmolDock::Wrapper {

    enum class ScoringFunctionEnum {
        VinaLike
    };

    std::function<double(const sd::iConformer &, const sd::iTransform &, const sd::iProtein &)>
    getScoringFunc(ScoringFunctionEnum which) {
        if (which == ScoringFunctionEnum::VinaLike) {
            return sd::Score::vina_like_rigid_inter_scoring_func;
        }

        return sd::Score::vina_like_rigid_inter_scoring_func; // By default return Vina like
    }

}


void export_ScoringFunctions() {
    p::class_<std::function<double(const sd::iConformer &, const sd::iTransform &, const sd::iProtein &)> >(
            "ScoringFunction")
        //.def("populateFromPDB", &sd::Molecule::populateFromPDB)
        //.def("set", &World::set)
        //.add_property("rovalue", &Num::get)
        //.add_property("value", &Num::get, &Num::set);
            ;
    p::def("getScoringFunc", &sd::Wrapper::getScoringFunc);

    p::enum_<sd::Wrapper::ScoringFunctionEnum>("ScoringFuncType")
            .value("VinaLike", sd::Wrapper::ScoringFunctionEnum::VinaLike);


}