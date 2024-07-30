///////////////////////////////////////////////////////////////////////
//  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
///////////////////////////////////////////////////////////////////////


#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

#define UINT32      unsigned int
#define INT32       int
#define UINT64      unsigned long long
// #define COUNTER     unsigned long long


// #define NOT_TAKEN 0
// #define TAKEN 1

#define FAILURE 0
#define SUCCESS 1

//JD2_2_2016
//typedef enum {
//  OPTYPE_OP               =2,
//  OPTYPE_BRANCH_COND      =3,
//  OPTYPE_RET              =4,
//  OPTYPE_BRANCH           =6,
//  OPTYPE_INDIRECT         =7,
//  OPTYPE_MAX              =8
//}OpType;

//JD2_17_2016 break down types into COND/UNCOND
typedef enum {
  OPTYPE_OP               =2,

  OPTYPE_RET_UNCOND,
  OPTYPE_JMP_DIRECT_UNCOND,
  OPTYPE_JMP_INDIRECT_UNCOND,
  OPTYPE_CALL_DIRECT_UNCOND,
  OPTYPE_CALL_INDIRECT_UNCOND,

  OPTYPE_RET_COND,
  OPTYPE_JMP_DIRECT_COND,
  OPTYPE_JMP_INDIRECT_COND,
  OPTYPE_CALL_DIRECT_COND,
  OPTYPE_CALL_INDIRECT_COND,

  OPTYPE_ERROR,

  OPTYPE_MAX
}OpType;



typedef enum {
    NoBranch,
    Return,
    CallDirect,
    CallIndirect,
    DirectCond,
    DirectUncond,
    IndirectCond,
    IndirectUncond,
    MAX
} BrType;


/** Some helper functions */
OpType convertBrType(BrType type);

#define PRINTDEBUG 0

#define DPRINTF(...) \
    if (PRINTDEBUG) [[unlikely]] { \
        printf(__VA_ARGS__); \
    }

#define DPRINTIF(cond, ...) \
    if (PRINTDEBUG && (cond)) [[unlikely]] { \
        printf(__VA_ARGS__); \
    }

#define PRINTIF(cond, ...) \
    if (cond) [[unlikely]] { \
        printf(__VA_ARGS__); \
    }

#endif

