//
// Created by eliane on 16/03/19.
//

#ifndef SMOLDOCK_CSVREADER_H
#define SMOLDOCK_CSVREADER_H

#include <string>
#include <vector>
#include <map>

namespace SmolDock {


    class CSVReader {
    public:
        CSVReader(std::string filename_, std::string delimiter_, bool hasHeader_, std::string escapeChar_ = "\\", std::string quoteChar_ = "\"");

        std::vector<std::map<std::string,std::string>> getRowsAsMap();

    private:
        std::string filename;
        std::string delimiter;
        bool hasHeader;
        std::string escapeChar;
        std::string quoteChar;
    };


}


#endif //SMOLDOCK_CSVREADER_H
