//
// Created by eliane on 02/12/18.
//

#ifndef SMOLDOCK_TIMINGSLOG_H
#define SMOLDOCK_TIMINGSLOG_H

#include <chrono>

#ifdef SMOLDOCK_VERBOSE_DEBUG
#define record_timings(A) auto A = std::chrono::system_clock::now()
#endif

#ifndef SMOLDOCK_VERBOSE_DEBUG
#define record_timings(A) ;;
#endif

#endif //SMOLDOCK_TIMINGSLOG_H
