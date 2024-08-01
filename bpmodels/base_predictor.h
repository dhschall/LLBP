
#pragma once
// #ifndef __BASE_PREDICTOR__
// #define __BASE_PREDICTOR__

#include "utils/common.h"

class BasePredictor {
    static inline UINT32 SatIncrement(UINT32 x, UINT32 max) {
        if (x < max) return x + 1;
        return x;
    }

    static inline UINT32 SatDecrement(UINT32 x) {
        if (x > 0) return x - 1;
        return x;
    }

   public:
    BasePredictor() {};
    virtual ~BasePredictor() = default;

    virtual bool GetPrediction(uint64_t PC) = 0;
    virtual void FirstTimeUpdate(uint64_t PC, bool taken,
                                uint64_t branchTarget) {};
    virtual void UpdatePredictor(uint64_t PC, bool resolveDir,
                                 bool predDir, uint64_t branchTarget) = 0;

    virtual void TrackOtherInst(uint64_t PC, OpType opType, bool taken,
                                uint64_t branchTarget) = 0;

    virtual void PrintStat(double NUMINST) {};
    virtual void DumpTables(std::string filename) {};
    virtual void LoadTables(std::string filename) {};
    virtual void StartTracer(std::string filename) {};
    virtual void tick() {};
    virtual void resetStats() {};
    virtual void btbMiss() {};
    virtual void setState(bool warmup=false) {};
};


BasePredictor* CreateBP(std::string bp_name);

// #endif  //__BASE_PREDICTOR__
