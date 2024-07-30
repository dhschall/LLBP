
#include "common.h"

using namespace std;



OpType convertBrType(BrType type) {
    switch (type) {
        case NoBranch:
            return OPTYPE_OP;
        case Return:
            return OPTYPE_RET_UNCOND;
        case CallDirect:
            return OPTYPE_CALL_DIRECT_UNCOND;
        case CallIndirect:
            return OPTYPE_CALL_INDIRECT_UNCOND;
        case DirectCond:
            return OPTYPE_JMP_DIRECT_COND;
        case DirectUncond:
            return OPTYPE_JMP_DIRECT_UNCOND;
        case IndirectCond:
            return OPTYPE_JMP_INDIRECT_COND;
        case IndirectUncond:
            return OPTYPE_JMP_INDIRECT_UNCOND;
        default:
            return OPTYPE_ERROR;
    }
}
