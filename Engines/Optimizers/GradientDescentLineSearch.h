//
// Created by eliane on 23/12/18.
//

#ifndef SMOLDOCK_GRADIENTDESCENTLINESEARCH_H
#define SMOLDOCK_GRADIENTDESCENTLINESEARCH_H

#include <functional>

#include <Engines/Math/GradientCalculator.h>


#include "Engines/Internals/iConformer.h"
#include "Engines/Internals/iTransform.h"
#include "Engines/Internals/iProtein.h"


namespace SmolDock {

    class GradientDescentLineSearch {
    public:
        GradientDescentLineSearch(std::function<double(const iConformer&, const iTransform&,const iProtein&)>  scorFunc,
            double differentialUpsilon = 0.01);

        void setProtein(const iProtein* prot);
        void setStartingConformer(const iConformer* conformer);
        bool optimize();

        iConformer getFinalConformer();
        double getScore();
        unsigned int getIterationNumber();

    private:
        std::function<double(const iConformer&, const iTransform&,const iProtein&)> scoringFunction;
        const iConformer* startingConformer;
        iConformer currentConformer;
        const iProtein* protein;

        iTransform currentTransform;
        double score = 0;
        unsigned int iterationNbr = 0;

        double differential_epsilon;

        iConformer getCurrentConformerWithTransform();

    };

}


#endif //SMOLDOCK_GRADIENTDESCENTLINESEARCH_H
