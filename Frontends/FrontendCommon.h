//
// Created by briand on 3/21/19.
//

#ifndef SMOLDOCK_FRONTENDCOMMON_H
#define SMOLDOCK_FRONTENDCOMMON_H

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>

namespace SmolDock {

    inline void setupLogPrinting(bool noDebug = false, bool onlyError = false)
    {
        if(noDebug)
        {
            if(onlyError)
            {
                boost::log::core::get()->set_filter
                        (boost::log::trivial::severity >= boost::log::trivial::error);
            }else {
                boost::log::core::get()->set_filter
                        (boost::log::trivial::severity >= boost::log::trivial::info);
            }

        }else {
            boost::log::core::get()->set_filter
                    (boost::log::trivial::severity >= boost::log::trivial::debug);
        }


        auto console_logger = boost::log::add_console_log(std::cout);
        console_logger->set_formatter([](boost::log::record_view const &rec, boost::log::formatting_ostream &strm) {
            if (rec[boost::log::trivial::severity] == boost::log::trivial::trace) {
                strm << " T  "; //         use TRACE_LOG(); macro for auto file:line:function
            } else if (rec[boost::log::trivial::severity] == boost::log::trivial::debug) {
                strm << "{D} ";
            } else if (rec[boost::log::trivial::severity] == boost::log::trivial::info) {
                strm << "    ";
            } else if (rec[boost::log::trivial::severity] == boost::log::trivial::warning) {
                strm << "[!] ";
            } else if (rec[boost::log::trivial::severity] >= boost::log::trivial::error) {
                strm << "[E] ";
            }

            strm << rec[boost::log::expressions::smessage];
        });
    }

}

#endif //SMOLDOCK_FRONTENDCOMMON_H
