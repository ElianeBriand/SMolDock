

#include "Engines/Internals/iTransform.h"
#include "Engines/Internals/InternalsUtilityFunctions.h"

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include <Utilities/IntermediateConformerCollector.h>

SmolDock::IntermediateConformerCollector* conformerCollector;

#include <random>

#include <Eigen/Eigen>
#include <Vc/Vc>
#include <Engines/ScoringFunctions/VinaLikeCommon.h>

namespace sd = SmolDock;


BOOST_AUTO_TEST_SUITE(Engines_ts)

    BOOST_AUTO_TEST_SUITE(geometry_ts)

    BOOST_AUTO_TEST_SUITE_END();

    BOOST_AUTO_TEST_SUITE(scoreComponents_ts)


        struct ExampleDistanceFixture {
            ExampleDistanceFixture() :
            mt(452),
            dist(0.0, 12.0){
                //static_assert(Vc::Vector<double>::Size >= 3);

                for(unsigned int i = 0; i < numTestPoint; i++) {
                    double random_distance = dist(mt);
                    distances_array[i] = random_distance;
                    distances_memoryblock.scalar(i) = random_distance;
                    BOOST_TEST_MESSAGE("ExampleDistanceFixture distance[" + boost::lexical_cast<std::string>(i) + "] = " + boost::lexical_cast<std::string>(random_distance));
                }
            }

            ~ExampleDistanceFixture() = default;

            std::mt19937 mt;
            std::uniform_real_distribution<double> dist;

            static constexpr const unsigned int numTestPoint = Vc::Vector<double>::Size * 4; // Multiple of vect size makes code easier
            std::array<double, numTestPoint> distances_array;
            Vc::Memory<Vc::Vector<double>,numTestPoint> distances_memoryblock;
        };

        BOOST_FIXTURE_TEST_CASE(gaussVectorized_equivalency, ExampleDistanceFixture,
                                *boost::unit_test::tolerance(0.001)) {

            std::array<double, numTestPoint> result_nonvectorized_gauss1;
            std::array<double, numTestPoint> result_nonvectorized_gauss2;
            for(unsigned int i = 0; i < numTestPoint; i++) {
                result_nonvectorized_gauss1[i] = sd::Score::vinaGaussComponent(distances_array[i], 0.0, 0.5);
                result_nonvectorized_gauss2[i] = sd::Score::vinaGaussComponent(distances_array[i], 3.0, 2.0);
            }

            Vc::Memory<Vc::Vector<double>,numTestPoint> result_vectorized_gauss1;
            Vc::Memory<Vc::Vector<double>,numTestPoint> result_vectorized_gauss2;
            for (unsigned int i = 0; i < result_vectorized_gauss1.vectorsCount(); ++i) {
                Vc::Vector<double> distance_simdvector = distances_memoryblock.vector(i);

                result_vectorized_gauss1.vector(i) = sd::Score::vinaGaussComponent(distance_simdvector, 0.0, 0.5);
                result_vectorized_gauss2.vector(i) = sd::Score::vinaGaussComponent(distance_simdvector, 3.0, 2.0);
            }

            for(unsigned int i = 0; i < numTestPoint; i++) {
                BOOST_TEST(result_nonvectorized_gauss1[i] == result_vectorized_gauss1.scalar(i));
            }

            for(unsigned int i = 0; i < numTestPoint; i++) {
                BOOST_TEST(result_nonvectorized_gauss2[i] == result_vectorized_gauss2.scalar(i));
            }

        }

        BOOST_FIXTURE_TEST_CASE(repulsionVectorized_equivalency, ExampleDistanceFixture,
                                *boost::unit_test::tolerance(0.001)) {

            std::array<double, numTestPoint> result_nonvectorized_repulsion;
            for(unsigned int i = 0; i < numTestPoint; i++) {
                result_nonvectorized_repulsion[i] = sd::Score::vinaRepulsionComponent(distances_array[i], 0.0);
            }

            Vc::Memory<Vc::Vector<double>,numTestPoint> result_vectorized_repulsion;
            for (unsigned int i = 0; i < result_vectorized_repulsion.vectorsCount(); ++i) {
                Vc::Vector<double> distance_simdvector = distances_memoryblock.vector(i);

                result_vectorized_repulsion.vector(i) = sd::Score::vinaRepulsionComponent(distance_simdvector, 0.0);
            }

            for(unsigned int i = 0; i < numTestPoint; i++) {
                BOOST_TEST(result_nonvectorized_repulsion[i] == result_vectorized_repulsion.scalar(i));
            }


        }

    BOOST_AUTO_TEST_SUITE_END();


    BOOST_AUTO_TEST_CASE(internals_TypeProperties) {
        BOOST_TEST(std::is_polymorphic<SmolDock::iConformer>::value == false);

    }


BOOST_AUTO_TEST_SUITE_END();