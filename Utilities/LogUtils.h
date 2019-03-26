//
// Created by briand on 3/25/19.
//

#ifndef SMOLDOCK_LOGUTILS_H
#define SMOLDOCK_LOGUTILS_H

#include <iostream>
#include <sstream>

template<typename T>
inline std::string vectorToString(const std::vector<T> vector) {
    std::stringstream ss;

    ss << "[ ";
    if (vector.size() != 0) {
        for (unsigned int j = 0; j < vector.size() - 1; ++j) {
            ss << " " << vector[j] << " ,";
        }
        ss << " " << vector[vector.size() - 1];
    }
    ss << " ]";

    return ss.str();
}

#endif //SMOLDOCK_LOGUTILS_H
