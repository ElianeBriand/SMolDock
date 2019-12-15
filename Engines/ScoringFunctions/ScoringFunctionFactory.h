//
// Created by eliane on 08/01/19.
//

#ifndef SMOLDOCK_SCORINGFUNCTIONFACTORY_H
#define SMOLDOCK_SCORINGFUNCTIONFACTORY_H

#include <memory>

#include <Engines/Internals/iProtein.h>
#include <boost/serialization/serialization.hpp>

#include "ScoringFunctionInterface.h"

namespace SmolDock::Score {


    enum class ScoringFunctionType {
        VinaRigid,
        Vina,
        VinaCovalentReversible,

    };

    std::shared_ptr<ScoringFunction> scoringFunctionFactory(ScoringFunctionType t, const iConformer &conformer,
                                                            const iProtein &protein, const iTransform &transform,
                                                            double differential_upsilon,
                                                            bool useNonDefaultCoefficients = false);

    std::string scoringFunctionTypeToString(ScoringFunctionType t);

    inline std::ostream& operator<<(std::ostream & os, ScoringFunctionType & scoringFunctionType) {
         switch (scoringFunctionType) {
          case ScoringFunctionType::VinaRigid:
              os << "VinaRigid";
              break;
          case ScoringFunctionType::Vina:
              os << "Vina";
              break;
          case ScoringFunctionType::VinaCovalentReversible:
              os << "VinaCovalentReversible";
      }
      return os;
    }


}

namespace boost {
namespace serialization {

    template<class Archive>
    void serialize(Archive & ar, SmolDock::Score::ScoringFunctionType & type, const unsigned int version)
    {
        long unsigned int code;
        ar & code;
        type = ( SmolDock::Score::ScoringFunctionType) code;

    }

} // namespace serialization
} // namespace boost

#endif //SMOLDOCK_SCORINGFUNCTIONFACTORY_H


