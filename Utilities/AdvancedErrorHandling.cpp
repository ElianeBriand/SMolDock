//
// Created by briand on 3/28/19.
//

#include <iostream>


#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <cxxabi.h>
#include <iomanip>
#include <exception>
#include <typeinfo>
#include <type_traits>

#define UNW_LOCAL_ONLY

#include <libunwind.h>

namespace {


    static std::terminate_handler previousHandler;

    char* get_demangled_name(char const* const symbol) noexcept {
        if (!symbol) {
            char* strPtr = (char*) malloc(sizeof(char) * strlen("<null>") + 1);
            strcpy(strPtr, "<null>");
            return strPtr;
        }
        int status = -4;
        char* demangledName = abi::__cxa_demangle(symbol, nullptr, nullptr, &status);
        if (!demangledName) {
            char* strPtr = (char*) malloc(sizeof(char) * strlen(symbol) + 1 + sizeof(char) * strlen("?: ") + 1);
            strcpy(strPtr, "?: ");
            strcat(strPtr, symbol);
            return strPtr;
        }
        return demangledName;
    }

    void
    print_reg(std::ostream& _out, unw_word_t reg) noexcept {
        constexpr std::size_t address_width = std::numeric_limits<std::uintptr_t>::digits / 4;
        _out << "0x" << std::setfill('0') << std::setw(address_width) << reg;
    }

    char symbol[1024];


    void print_stacktrace(std::ostream& _out) noexcept {
        unw_cursor_t cursor;
        unw_context_t context;
        unw_getcontext(&context);
        unw_init_local(&cursor, &context);
        _out << std::hex << std::uppercase;
        while (0 < unw_step(&cursor)) {
            _out << "    ";
            unw_word_t ip = 0;
            unw_get_reg(&cursor, UNW_REG_IP, &ip);
            if (ip == 0) {
                break;
            }
            unw_word_t sp = 0;
            unw_get_reg(&cursor, UNW_REG_SP, &sp);
            print_reg(_out, ip);
            _out << ": (SP:";
            print_reg(_out, sp);
            _out << ") ";
            unw_word_t offset = 0;
            if (unw_get_proc_name(&cursor, symbol, sizeof(symbol), &offset) == 0) {
                char* demangledName = get_demangled_name(symbol);
                _out << "(" << demangledName << " + 0x" << offset << ")\n\n";
                free(demangledName);
            } else {
                _out << "-- error: unable to obtain symbol name for this frame\n\n";
            }
        }
        _out << std::flush;
    }

/* Obtain a backtrace and print it to stdout. */
    void print_trace() noexcept {
        print_stacktrace(std::cout);
    }


    [[noreturn]] void custom_terminate_handler() {

        static bool terminating;
        if (terminating) {
            std::cout << "Terminate handler called recursively\n" << std::endl;
            abort();
        }
        terminating = true;

        std::cout << "\n\n\n  ------------- TERMINATE HANDLER -------------" << std::endl;

        if (auto exc = std::current_exception()) {


            std::type_info* t = abi::__cxa_current_exception_type();
            char const* name = t->name();
            char* demangledName = get_demangled_name(name);

            std::cout << "  Exception : Unhandled exception of type " << demangledName << std::endl;


            free(demangledName);
        } else {
            std::cout << "  Exception : terminate called without active exception." << std::endl;
        }
        std::cout << "\n  Stack trace : " << std::endl;
        print_trace();
        std::cout << " ---------------------------------------------\n\n\n" << std::endl;
        std::abort();
    }


}

void setupAdvancedErrorHandling() {
    previousHandler = std::set_terminate(custom_terminate_handler);
}

namespace boost {
    void assertion_failed(char const* expr, char const* function,
                          char const* file, long line) {

        std::cout << "\n\n\n ------------- ASSERTION FAILED -------------" << std::endl;
        std::cout << "  Assertion : " << expr << std::endl;
        std::cout << "  Function  : " << function << std::endl;
        std::cout << "  File      : " << file << " @ line " << line << std::endl;
        std::cout << "\n  Stack trace : " << std::endl;
        print_trace();

        std::cout << " --------------------------------------------\n\n\n" << std::endl;

        std::abort();
    }


    void assertion_failed_msg(char const* expr, char const* msg,
                              char const* function, char const* file, long line) {


    }
}