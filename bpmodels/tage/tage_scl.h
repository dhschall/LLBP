/* MIT License
 *
 * Copyright (c) 2024 David Schall and EASE lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * The code is based on the 64KiB TAGE-SC-L branch predictor by Andre Seznec
 * provided in the CBP-5 competition.
 * It was reformated and made easier to configure. Furthermore, the code
 * adds a lot of statistics and debugging information to the predictor.
 */

#pragma once

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "tage.h"

//----------------------
// #define TAGEInf
// #define TAGE8k
//----------------------

namespace LLBP {


struct TSCLConfig;

class TageSCL : public TageBase {
   public:
    TageSCL(TSCLConfig config);
    ~TageSCL();

    void DumpTables(std::string filename) override;
    void PrintStat(double NUMINST) override;

   protected:
    typedef TageBase Tage;

    // Override some base class functions
    bool predict(uint64_t pc) override;
    void updateTables(uint64_t pc, bool resolveDir, bool predDir) override;

    void updateStats(bool taken, bool predtaken, uint64_t PC) override;
    void updateHistory(const uint64_t pc, const bool taken, const OpType opType,
                       const uint64_t branchTarget) override;

    void init_predictor();
    void predictorsize();
    void resetStats() override {
        Tage::resetStats();
        sclstats = {};
    };
    // Internal prediction of the SCL part
    bool SCLpredict(uint64_t pc);




// --- SC --- //
// The max number is not the actual number of entries just for the
// initialization

#ifdef TAGEInf


// Bias tables
#define LogBiasMax 20
    int8_t Bias[(1 << LogBiasMax)];
    int8_t BiasSK[(1 << LogBiasMax)];
    int8_t BiasBank[(1 << LogBiasMax)];

    // In all th GEHL components, the two tables with the shortest history
    // lengths have only half of the entries.

    // IMLI-SIC -> Micro 2015  paper: a big disappointment on  CBP2016 traces
    long long IMLIcount;  // use to monitor the iteration number

// Not used by INF SCL ------------------------------
#define LOGBWNB 7
#define BWNB 2
    int BWm[BWNB] = { 16, 8 };
    int8_t BWGEHLA[BWNB][(1 << LOGBWNB)] = { {0} };
    int8_t *BWGEHL[BWNB];
    long long BWHIST;


#define LOGINB 19  // 128-entry
#define INB 1
    int Im[INB] = {8};
    int8_t IGEHLA[INB][(1 << LOGINB)] = {{0}};
    int8_t *IGEHL[INB];

// global branch GEHL
#define LOGGNB 19  // 1 1K + 2 * 512-entry tables
#define GNB 3
    int Gm[GNB] = {40, 24, 10};
    int8_t GGEHLA[GNB][(1 << LOGGNB)] = {{0}};
    int8_t *GGEHL[GNB];

// first local history
#define LOGLNB 19  // 1 1K + 2 * 512-entry tables
#define LNB 3
    int Lm[LNB] = {11, 6, 3};
    int8_t LGEHLA[LNB][(1 << LOGLNB)] = {{0}};
    int8_t *LGEHL[LNB];
#define LOGLOCAL 19
#define NLOCAL (1 << LOGLOCAL)
#define INDLOCAL ((PC ^ (PC >> 2)) & (NLOCAL - 1))
    long long L_shist[NLOCAL];  // local histories


// Variation on the IMLI predictor
#define LOGIMNB 19  // 2* 256 -entry
#define IMNB 2
    int IMm[IMNB] = {10, 4};
    int8_t IMGEHLA[IMNB][(1 << LOGIMNB)] = {{0}};
    int8_t *IMGEHL[IMNB];
    long long IMHIST[256];


// variation on global branch history
#define PNB 3
#define LOGPNB 19  // 1 1K + 2 * 512-entry tables
    int Pm[PNB] = {25, 16, 9};
    int8_t PGEHLA[PNB][(1 << LOGPNB)] = {{0}};
    int8_t *PGEHL[PNB];


// second local history
#define LOGSNB 19  // 1 1K + 2 * 512-entry tables
#define SNB 3
    int Sm[SNB] = {16, 11, 6};
    int8_t SGEHLA[SNB][(1 << LOGSNB)] = {{0}};

    int8_t *SGEHL[SNB];
#define LOGSECLOCAL 4
#define NSECLOCAL (1 << LOGSECLOCAL)  // Number of second local histories
#define INDSLOCAL (((PC ^ (PC >> 5))) & (NSECLOCAL - 1))
    long long S_slhist[NSECLOCAL];

// third local history
#define LOGTNB 19  // 2 * 512-entry tables
#define TNB 2
    int Tm[TNB] = {9, 4};
    int8_t TGEHLA[TNB][(1 << LOGTNB)] = {{0}};

    int8_t *TGEHL[TNB];
#define NTLOCAL 19
#define INDTLOCAL                \
    (((PC ^ (PC >> (LOGTNB)))) & \
     (NTLOCAL - 1))  // different hash for the history
    long long T_slhist[NTLOCAL];


#else // No TAGEInf ------------------------------



// Bias tables
#define LogBiasMax 15
    int8_t Bias[(1 << LogBiasMax)];
    int8_t BiasSK[(1 << LogBiasMax)];
    int8_t BiasBank[(1 << LogBiasMax)];

    // In all th GEHL components, the two tables with the shortest history
    // lengths have only half of the entries.

    // IMLI-SIC -> Micro 2015  paper: a big disappointment on  CBP2016 traces
    long long IMLIcount;  // use to monitor the iteration number


#ifdef TAGE8k

#define LOGINB 7
#define INB 1
    int Im[INB] = { 8 };
    int8_t IGEHLA[INB][(1 << LOGINB)] = { {0} };
    int8_t *IGEHL[INB];

//global branch GEHL
#define LOGGNB 7
#define GNB 2
    int Gm[GNB] = {6,3};
    int8_t GGEHLA[GNB][(1 << LOGGNB)] = { {0} };
    int8_t *GGEHL[GNB];


//large local history
#define LOGLNB 7
#define LNB 2
    int Lm[LNB] = { 6, 3 };
    int8_t LGEHLA[LNB][(1 << LOGLNB)] = { {0} };
    int8_t *LGEHL[LNB];

#define  LOGLOCAL 6
#define NLOCAL (1<<LOGLOCAL)
#define INDLOCAL ((PC ^ (PC >>2)) & (NLOCAL-1))
    long long L_shist[NLOCAL];


#else // TAGE64k

#define LOGINB 8  // 128-entry
#define INB 1
    int Im[INB] = {8};
    int8_t IGEHLA[INB][(1 << LOGINB)] = {{0}};
    int8_t *IGEHL[INB];

// global branch GEHL
#define LOGGNB 10  // 1 1K + 2 * 512-entry tables
#define GNB 3
    int Gm[GNB] = {40, 24, 10};
    int8_t GGEHLA[GNB][(1 << LOGGNB)] = {{0}};
    int8_t *GGEHL[GNB];

// first local history
#define LOGLNB 10  // 1 1K + 2 * 512-entry tables
#define LNB 3
    int Lm[LNB] = {11, 6, 3};
    int8_t LGEHLA[LNB][(1 << LOGLNB)] = {{0}};
    int8_t *LGEHL[LNB];
#define LOGLOCAL 8
#define NLOCAL (1 << LOGLOCAL)
#define INDLOCAL ((PC ^ (PC >> 2)) & (NLOCAL - 1))
    long long L_shist[NLOCAL];  // local histories

#endif


// Only used by 8k TAGE ------------------------------
// Backward branch history
#define LOGBWNB 7
#define BWNB 2
    int BWm[BWNB] = { 16, 8 };
    int8_t BWGEHLA[BWNB][(1 << LOGBWNB)] = { {0} };
    int8_t *BWGEHL[BWNB];
    long long BWHIST;



// Only used by 64k TAGE ------------------------------
// Variation on the IMLI predictor
#define LOGIMNB 9  // 2* 256 -entry
#define IMNB 2
    int IMm[IMNB] = {10, 4};
    int8_t IMGEHLA[IMNB][(1 << LOGIMNB)] = {{0}};
    int8_t *IMGEHL[IMNB];
    long long IMHIST[256];


// variation on global branch history
#define PNB 3
#define LOGPNB 9  // 1 1K + 2 * 512-entry tables
    int Pm[PNB] = {25, 16, 9};
    int8_t PGEHLA[PNB][(1 << LOGPNB)] = {{0}};
    int8_t *PGEHL[PNB];


// second local history
#define LOGSNB 9  // 1 1K + 2 * 512-entry tables
#define SNB 3
    int Sm[SNB] = {16, 11, 6};
    int8_t SGEHLA[SNB][(1 << LOGSNB)] = {{0}};

    int8_t *SGEHL[SNB];
#define LOGSECLOCAL 4
#define NSECLOCAL (1 << LOGSECLOCAL)  // Number of second local histories
#define INDSLOCAL (((PC ^ (PC >> 5))) & (NSECLOCAL - 1))
    long long S_slhist[NSECLOCAL];

// third local history
#define LOGTNB 10  // 2 * 512-entry tables
#define TNB 2
    int Tm[TNB] = {9, 4};
    int8_t TGEHLA[TNB][(1 << LOGTNB)] = {{0}};

    int8_t *TGEHL[TNB];
#define NTLOCAL 16
#define INDTLOCAL                \
    (((PC ^ (PC >> (LOGTNB)))) & \
     (NTLOCAL - 1))  // different hash for the history
    long long T_slhist[NTLOCAL];


#endif

#define LOGSIZEUPMAX 6  // not worth increasing
    int updatethreshold;
    int Pupdatethreshold[(1 << LOGSIZEUPMAX)];  // size is fixed by LOGSIZEUP
    int8_t WG[(1 << LOGSIZEUPMAX)];
    int8_t WL[(1 << LOGSIZEUPMAX)];
    int8_t WS[(1 << LOGSIZEUPMAX)];
    int8_t WT[(1 << LOGSIZEUPMAX)];
    int8_t WP[(1 << LOGSIZEUPMAX)];
    int8_t WI[(1 << LOGSIZEUPMAX)];
    int8_t WIM[(1 << LOGSIZEUPMAX)];
    int8_t WB[(1 << LOGSIZEUPMAX)];
    int8_t WBW[(1 << LOGSIZEUPMAX)];

    long long GHIST;
    int LSUM;
    bool sc_pred, tage_scl_pred;

    // The thereshold for the SC prediction
    int THRES;
    //
    int64_t gIndex(uint64_t PC, int64_t bhist, int logs, int nbr, int i);

    int Gpredict(uint64_t PC, long long BHIST, int *length, int8_t **tab,
                 int NBR, int logs, int8_t *W);
    void Gupdate(uint64_t PC, bool taken, long long BHIST, int *length,
                 int8_t **tab, int NBR, int logs, int8_t *W);

    void SCpredict(uint64_t pc);
    void SCUpdate(uint64_t PC, bool resolveDir, bool predDir);

    void SCLPredict(uint64_t pc);
    void SCLUpdate(uint64_t PC, bool resolveDir, bool predDir);


    // Two hooks to add additional information into the SC prediction
    virtual int compPartial(uint64_t pc) { return 0; }
    virtual void updatePartial(uint64_t PC, bool resolveDir) {}

    // ---- LOOP PREDICTOR --- //
    class lentry  // loop predictor entry
    {
       public:
        uint16_t NbIter;       // 10 bits
        uint8_t confid;        // 4bits
        uint16_t CurrentIter;  // 10 bits

        uint16_t TAG;  // 10 bits
        uint8_t age;   // 4 bits
        bool dir;      // 1 bit

        // 39 bits per entry
        lentry() {
            confid = 0;
            CurrentIter = 0;
            NbIter = 0;
            TAG = 0;
            age = 0;
            dir = false;
        }
    };
    const int LogL;  // log of number of entries in the loop predictor
    lentry *ltable;  // loop predictor table
    // variables for the loop predictor
    bool loop_pred;  // loop predictor prediction
    int LIB;
    int LI;
    int LHIT;         // hitting way in the loop predictor
    int LTAG;         // tag on the loop predictor
    bool LVALID;      // validity of the loop predictor prediction
    int8_t WITHLOOP;  // counter to monitor whether or not loop prediction is
                      // beneficial

    int lindex(uint64_t PC);
    bool getloop(uint64_t PC);
    void loopupdate(uint64_t PC, bool Taken, bool ALLOC);
    int _lSeed;
    // just a simple pseudo random number generator: use available information
    //  to allocate entries  in the loop predictor
    int MYLoopRANDOM() {
        _lSeed++;
        // Seed ^= phist;
        _lSeed = (_lSeed >> 21) + (_lSeed << 11);
        _lSeed = (_lSeed >> 10) + (_lSeed << 22);
        return (_lSeed);
    };

    // Chooser ---------------------
    enum {
      LOOP = Tage::LAST_TAGE_PROVIDER_TYPE +1,
      STC = LOOP + 1,
      LAST_SCL_PROVIDER_TYPE = STC
    };
    int scl_provider;

    // The two counters used to choose between TAGE and SC on Low Conf SC
    int8_t FirstH, SecondH;
    const int ConfWidth = 7;  // for the counters in the choser

    const bool useSC;
    const bool useLoop;

    // Disable the chooser confidence counter
    const bool disableConfCounter;

    // Some stats
    struct
    {
        /* data */
        int provLoop = 0;
        int loopCorrect = 0;
        int loopIncorrect = 0;
        int provSC = 0;
        int scCorrect = 0;
        int scIncorrect = 0;
        int provTage = 0;
        int tageCorrect = 0;
        int tageIncorrect = 0;
        int tageMisses = 0;
        int tsclMisses = 0;
    } sclstats;

};


// Configurations
struct TSCLConfig {
    TageConfig tageConfig;
    bool useSC = true;
    bool useLoop = true;
    int LogL = 5;

    void print() const {
        printf("SCL Config: useSC=%d, useLoop=%d\n", useSC, useLoop);
    }
};
inline const TSCLConfig TSCL64kCfgDefault = {};


// Configurations
class TageSCL64k : public TageSCL {
   public:
    TageSCL64k(void)
        : TageSCL(TSCL64kCfgDefault)
    {}
};

class TageSCL512k : public TageSCL {
   public:
    TageSCL512k(void)
        : TageSCL(TSCLConfig
            {
                .tageConfig = TageConfig  {
                    .LogG = 13,    // 64K (2**10 entries) * 8
                }
            })
    {}
};


};  // namespace LLBP
