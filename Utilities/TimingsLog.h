/*
 * Copyright (c) 2018 Eliane Briand
 *
 * This file is part of SmolDock.
 *
 * SmolDock is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef SMOLDOCK_TIMINGSLOG_H
#define SMOLDOCK_TIMINGSLOG_H

#include <chrono>

#ifdef SMOLDOCK_VERBOSE_DEBUG
#define record_timings(A) auto A = std::chrono::system_clock::now()
#define TRACE_LOG()  BOOST_LOG_TRIVIAL(trace) << __FUNCTION__ << " <" << __FILE__ << ":" << __LINE__ << ":" <<"> "
#endif

#ifndef SMOLDOCK_VERBOSE_DEBUG
#define record_timings(A) ;;
#define TRACE_LOG()  ;;
#endif

#endif //SMOLDOCK_TIMINGSLOG_H