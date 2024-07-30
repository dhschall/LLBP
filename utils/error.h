


#include <cassert>
// #define assertm(exp, ...) assert(((void)msg, exp))

#define panic_if(cond, ...)                                  \
    do {                                                     \
        if (cond) {                                       \
            printf(__VA_ARGS__);                          \
            assert(0);                                      \
        }                                                    \
    } while (0)
