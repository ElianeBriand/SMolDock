//
// Created by eliane on 16/03/19.
//

#include "CSVReader.h"

#include <iostream>
#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>


namespace SmolDock {

    CSVReader::CSVReader(std::string filename_, std::string delimiter_, bool hasHeader_, std::string escapeChar_,
                         std::string quoteChar_) :
                         filename(filename_),
                         delimiter(delimiter_),
                         hasHeader(hasHeader_),
                         escapeChar(escapeChar_),
                         quoteChar(quoteChar_)
                         {
    }

    std::vector<std::map<std::string,std::string>> CSVReader::getRowsAsMap() {
        std::vector<std::map<std::string,std::string>> ret;

        std::ifstream file(this->filename);

        std::vector<std::string> headers;

        if(this->hasHeader)
        {
            std::string         firstline;
            std::getline(file, firstline);

            boost::tokenizer<boost::escaped_list_separator<char> > tk(
                    firstline, boost::escaped_list_separator<char>('\\', '\t', '\"'));
            for (boost::tokenizer<boost::escaped_list_separator<char> >::iterator i(tk.begin());
                    i!=tk.end();++i)
            {
                headers.push_back(*i);
            }
        }else {
            std::ifstream file2(this->filename);
            std::string         firstline;
            std::getline(file2, firstline);

            unsigned int numColumn = 0;

            boost::tokenizer<boost::escaped_list_separator<char> > tk(
                    firstline, boost::escaped_list_separator<char>('\\', '\t', '\"'));
            for (boost::tokenizer<boost::escaped_list_separator<char> >::iterator i(tk.begin());
                 i!=tk.end();++i)
            {
                headers.push_back(boost::lexical_cast<std::string>(numColumn));
                numColumn++;
            }
        }

        unsigned int lineCount = 1;
        while(file)
        {
            std::string         line;
            std::getline(file, line);

            std::vector<std::string> row;
            boost::tokenizer<boost::escaped_list_separator<char> > tk(
                    line, boost::escaped_list_separator<char>('\\', '\t', '\"'));
            for (boost::tokenizer<boost::escaped_list_separator<char> >::iterator i(tk.begin());
                 i!=tk.end();++i)
            {
                row.push_back(*i);
            }

            if(row.size() == 0)
                continue;

            if(row.size() != headers.size())
            {
                BOOST_LOG_TRIVIAL(error) << "Encountered a row with " << row.size() << " column but header/first line had " << headers.size() << " columns.";
                BOOST_LOG_TRIVIAL(error) << "Line number " << lineCount <<" (not counting header) will thus be discarded";
                continue;
            }

            std::map<std::string,std::string> rowMap;
            for (unsigned int j = 0; j < row.size(); ++j) {
                rowMap[ headers.at(j) ] = row[j];
            }
            ret.push_back(rowMap);

            lineCount++;
        }

        BOOST_LOG_TRIVIAL(info) << "Read " << lineCount << " lines from " << this->filename;

        return ret;
    }


}