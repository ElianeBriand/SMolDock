

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

        BOOST_AUTO_TEST_CASE(transformApplicationFunction, *boost::unit_test::tolerance(0.001)) {

            Eigen::Vector3d Vector1 = {14.2, -4.5, 3};
            Eigen::Vector3d Vector2 = {-5.3, 4.8, -2.7};
            Eigen::Vector3d Vector3 = {0.0, 8.4, 7.3};

            Vc::Memory<Vc::Vector<double>, 3> x, y, z;

            x.scalar(0) = Vector1.x();
            y.scalar(0) = Vector1.y();
            z.scalar(0) = Vector1.z();

            x.scalar(1) = Vector2.x();
            y.scalar(1) = Vector2.y();
            z.scalar(1) = Vector2.z();

            x.scalar(2) = Vector3.x();
            y.scalar(2) = Vector3.y();
            z.scalar(2) = Vector3.z();

            std::cout << "x_init =  " << x << std::endl;
            std::cout << "y_init =  " << y << std::endl;
            std::cout << "z_init =  " << z << std::endl;

            std::cout << "v1_init = " << Vector1 << std::endl;
            std::cout << "v2_init = " << Vector2 << std::endl;
            std::cout << "v3_init = " << Vector3 << std::endl;

            Eigen::Vector3d transl;
            transl.x() = 2;
            transl.y() = -3;
            transl.z() = 6;


            Eigen::Quaternion<double> rota;
            rota.x() = 2;
            rota.y() = 3;
            rota.z() = -1;
            rota.w() = 2.3;

            SmolDock::iTransform tr;
            tr.rota = rota;
            tr.transl = transl;
            tr.doHousekeeping();

            rota.normalize();

            Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotMatrix = rota.toRotationMatrix();


            Eigen::Vector3d Res1 = (rotMatrix * Vector1) + transl;

            Eigen::Vector3d Res2 = (rotMatrix * Vector2) + transl;

            Eigen::Vector3d Res3 = (rotMatrix * Vector3) + transl;

            Eigen::Vector3d Res1_function = Vector1;
            Eigen::Vector3d Res2_function = Vector2;
            Eigen::Vector3d Res3_function = Vector3;

            applyRigidTransformInPlace(Res1_function, tr);
            applyRigidTransformInPlace(Res2_function, tr);
            applyRigidTransformInPlace(Res3_function, tr);


            const Vc::Vector<double> xold = x.vector(0);
            const Vc::Vector<double> yold = y.vector(0);
            const Vc::Vector<double> zold = z.vector(0);


            /*
             * rotMatrix =
             * | a1 a2 a3 |
             * | b1 b2 b3 |
             * | c1 c2 c3 |
             */
            const double a1 = rotMatrix.data()[0];
            const double a2 = rotMatrix.data()[1];
            const double a3 = rotMatrix.data()[2];
            const double b1 = rotMatrix.data()[3];
            const double b2 = rotMatrix.data()[4];
            const double b3 = rotMatrix.data()[5];
            const double c1 = rotMatrix.data()[6];
            const double c2 = rotMatrix.data()[7];
            const double c3 = rotMatrix.data()[8];

            x.vector(0) = a1 * xold + a2 * yold + a3 * zold;
            y.vector(0) = b1 * xold + b2 * yold + b3 * zold;
            z.vector(0) = c1 * xold + c2 * yold + c3 * zold;

            x = x.vector(0) + Vc::Vector<double>(Vc::One) * transl.x();
            y = y.vector(0) + Vc::Vector<double>(Vc::One) * transl.y();
            z = z.vector(0) + Vc::Vector<double>(Vc::One) * transl.z();

            std::cout << "x_res  =  " << x << std::endl;
            std::cout << "y_res  =  " << y << std::endl;
            std::cout << "z_res  =  " << z << std::endl;

            std::cout << "v1_res =  " << Res1 << std::endl;
            std::cout << "v2_res =  " << Res2 << std::endl;
            std::cout << "v3_res =  " << Res3 << std::endl;

            Vc::Vector<double> x_func = xold;
            Vc::Vector<double> y_func = yold;
            Vc::Vector<double> z_func = zold;

            applyRigidTransformInPlace(x_func,
                                       y_func,
                                       z_func,
                                       tr);

            std::cout << " == Manual matrix+trans , comparing Eigen::v3d and Vc::Vector" << std::endl;

            std::cout << "x_res  =  " << x << std::endl;
            std::cout << "y_res  =  " << y << std::endl;
            std::cout << "z_res  =  " << z << std::endl;
            std::cout << "v1_res =  " << Res1 << std::endl;
            std::cout << "v2_res =  " << Res2 << std::endl;
            std::cout << "v3_res =  " << Res3 << std::endl;

            BOOST_TEST(x.scalar(0) == Res1.x());
            BOOST_TEST(y.scalar(0) == Res1.y());
            BOOST_TEST(z.scalar(0) == Res1.z());

            BOOST_TEST(x.scalar(1) == Res2.x());
            BOOST_TEST(y.scalar(1) == Res2.y());
            BOOST_TEST(z.scalar(1) == Res2.z());

            BOOST_TEST(x.scalar(2) == Res3.x());
            BOOST_TEST(y.scalar(2) == Res3.y());
            BOOST_TEST(z.scalar(2) == Res3.z());


            std::cout << " == Manual matrix+trans vs applyRigidTransformInPlace, on Vc::Vector ==" << std::endl;

            std::cout << "x_res  =  " << x << std::endl;
            std::cout << "y_res  =  " << y << std::endl;
            std::cout << "z_res  =  " << z << std::endl;

            std::cout << "x_func_res =  " << x_func << std::endl;
            std::cout << "y_func_res =  " << y_func << std::endl;
            std::cout << "z_func_res =  " << z_func << std::endl;

            BOOST_TEST(x.scalar(0) == x_func[0]);
            BOOST_TEST(y.scalar(0) == y_func[0]);
            BOOST_TEST(z.scalar(0) == z_func[0]);

            BOOST_TEST(x.scalar(1) == x_func[1]);
            BOOST_TEST(y.scalar(1) == y_func[1]);
            BOOST_TEST(z.scalar(1) == z_func[1]);

            BOOST_TEST(x.scalar(2) == x_func[2]);
            BOOST_TEST(y.scalar(2) == y_func[2]);
            BOOST_TEST(z.scalar(2) == z_func[2]);


            std::cout << " == Manual matrix+trans vs applyRigidTransformInPlace, on Eigen::v3d ==" << std::endl;
            std::cout << "v1_res =  " << Res1 << std::endl;
            std::cout << "v2_res =  " << Res2 << std::endl;
            std::cout << "v3_res =  " << Res3 << std::endl;

            std::cout << "v1_res_func =  " << Res1_function << std::endl;
            std::cout << "v2_res_func =  " << Res2_function << std::endl;
            std::cout << "v3_res_func =  " << Res3_function << std::endl;


            BOOST_TEST(Res1_function.x() == Res1.x());
            BOOST_TEST(Res1_function.y() == Res1.y());
            BOOST_TEST(Res1_function.z() == Res1.z());

            BOOST_TEST(Res2_function.x() == Res2.x());
            BOOST_TEST(Res2_function.y() == Res2.y());
            BOOST_TEST(Res2_function.z() == Res2.z());

            BOOST_TEST(Res3_function.x() == Res3.x());
            BOOST_TEST(Res3_function.y() == Res3.y());
            BOOST_TEST(Res3_function.z() == Res3.z());


            //BOOST_CHECK((resVect3.z - (3.07758)) < 0.001);


        }

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