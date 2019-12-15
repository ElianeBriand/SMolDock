/*
 * Copyright (c) 2019 Eliane Briand
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

#include <stdexcept>

#include <boost/log/trivial.hpp>

#include "ExtractProteinFromBox.hpp"

#include <Engines/Internals/iProtein.h>
#include <Engines/AbstractDockingEngine.h>

namespace SmolDock::Engine {

	iProtein extractIProteinFromBoxSetting(Protein* protein, const AbstractDockingEngine::DockingBoxSetting& setting) {
			iProtein ret_iprot;
            if (setting.shape == AbstractDockingEngine::DockingBoxSetting::Shape::sphere) {
                ret_iprot = protein->getPartialiProtein_sphere(setting.center,setting.radius, 2.0);
            } else if (setting.shape == AbstractDockingEngine::DockingBoxSetting::Shape::cube) {
            	// Unsupported.
            	BOOST_LOG_TRIVIAL(error) << "Cube docking box unsupported yet";
            	throw std::invalid_argument("Cube docking box unsupported yet");
        	} else if (setting.shape == AbstractDockingEngine::DockingBoxSetting::Shape::whole_protein) {
        		ret_iprot = protein->getiProtein();
            }else{
            	BOOST_LOG_TRIVIAL(warning) << "Unknown docking box type: defaulting to whole_protein";
                ret_iprot = protein->getiProtein();
            }
            return ret_iprot;
	}

}