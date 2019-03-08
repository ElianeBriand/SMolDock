//
// Created by eliane on 03/03/19.
//

#ifndef SMOLDOCK_VERSION_H
#define SMOLDOCK_VERSION_H

#include <string>

namespace SmolDock {

    class Version {
    public:

        constexpr static unsigned int major = 0;
        constexpr static unsigned int minor = 2;
        constexpr static unsigned int patch = 1;

        constexpr static char const* tag = "";


    };

    inline std::string getVersionString() {
        std::string s;
        s = "v" + std::to_string(Version::major) + "." +
            std::to_string(Version::minor) + "." +
            std::to_string(Version::patch);
        if (!std::string(Version::tag).empty()) {
            s += "-" + std::string(Version::tag);
        }
        return s;
    }

}
#endif //SMOLDOCK_VERSION_H
