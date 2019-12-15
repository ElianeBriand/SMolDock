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

#include <iostream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// Test dependencies
#include <UnitTests/Fixtures/TripeptideFixture.hpp>


// Module under test
#include <Engines/Utils/ExtractProteinFromBox.hpp>

namespace sd = SmolDock;


BOOST_AUTO_TEST_SUITE(Engines_ts)

    BOOST_AUTO_TEST_SUITE(Utilities_ts)

        BOOST_FIXTURE_TEST_CASE(DockingBoxUtils_selection_test, TripeptideFixture, *boost::unit_test::tolerance(0.001)) {

            sd::Engine::AbstractDockingEngine::DockingBoxSetting setting;
            setting.shape = sd::Engine::AbstractDockingEngine::DockingBoxSetting::Shape::sphere;
            setting.center = {7.269, 25.291, 22.254};
            setting.radius = 1.0;

            sd::iProtein res = sd::Engine::extractIProteinFromBoxSetting(tripeptide_Protein.get(), setting);
            BOOST_TEST(res.x.size() == 6);
            BOOST_TEST(res.y.size() == 6);
            BOOST_TEST(res.z.size() == 6);
            BOOST_TEST(res.atomicRadius.size() == 6);
            BOOST_TEST(res.type.size() == 6);
            BOOST_TEST(res.variant.size() == 6);

            setting.radius = 2.0;
            res = sd::Engine::extractIProteinFromBoxSetting(tripeptide_Protein.get(), setting);
            BOOST_TEST(res.x.size() == 22);
            BOOST_TEST(res.y.size() == 22);
            BOOST_TEST(res.z.size() == 22);
            BOOST_TEST(res.atomicRadius.size() == 22);
            BOOST_TEST(res.type.size() == 22);
            BOOST_TEST(res.variant.size() == 22);



        }

    BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();