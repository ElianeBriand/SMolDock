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

#include <Frontends/CLI/FileInputOutput.hpp>


#include <boost/filesystem.hpp>

namespace SmolDock {

	ReceptorFiletype receptorFiletypeFromPath(const std::string& filename) {
		boost::filesystem::path p(filename);
		const std::string extension = p.extension();

		if(extension == ".pdb")
			return ReceptorFiletype::pdb;
		else
			return ReceptorFiletype::unsupported;
	}

}