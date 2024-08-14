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
 * Credit:
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
#include <unordered_map>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stack>

#include "bpmodels/base_predictor.h"
#include "bpmodels/components/hist_registers.h"



namespace LLBP {



struct TageConfig;

class TageBase : public BasePredictor {
   public:
    TageBase(TageConfig config);
    ~TageBase();

    bool GetPrediction(uint64_t PC) override;
    void FirstTimeUpdate(uint64_t PC, bool taken,
                        uint64_t branchTarget) override;
    void UpdatePredictor(uint64_t PC, bool resolveDir,
                         bool predDir, uint64_t branchTarget) override;
    void TrackOtherInst(uint64_t PC, OpType opType, bool taken,
                        uint64_t branchTarget) override;

    virtual void PrintStat(double instr) override;
    void tick() override {
        ticks++;
    };
    void resetStats() override;



   protected:
    static const int histbufferlen =
        8192;  // we use a 4K entries history buffer to store the branch history
               // (this allows us to explore using history length up to 4K)

    static const int MAXNHIST = 40;  // Constant limit for the number of tables

    const int nhist;  // twice the number of different histories

    const int nbanklow;  // number of banks in the shared bank-interleaved for the low
             // history lengths
    const int nbankhigh;  // number of banks in the shared bank-interleaved for the  history
             // lengths
    const int born;  // below born in the table for low history lengths, >= born in the
             // table for high history lengths,
    // we use 2-way associativity for the medium history lengths
    const int assoc_start;  // 2 -way assoc for those banks 0.4 %
    const int assoc_end;

    int SizeTable[MAXNHIST];

    /*in practice 2 bits or 3 bits par branch: around 1200 cond. branchs*/

    const int minhist;  // 6  // not optimized so far
    const int maxhist;  // 3000

    const int hystshift = 2;  // bimodal hysteresis shared by 4 entries
    const int LogB;  // 13      // log of number of entries in bimodal predictor
    const int
        LogG;  // 10 /* logsize of the  banks in the  tagged TAGE tables */
    const int Tbits;  // 8
                      // minimum width of the tags  (low history lengths), +4
                      // for high history
    // lengths

    bool NOSKIP[MAXNHIST];  // to manage the associativity for different
                             // history lengths

    const int nnn =
        1;  // number of extra entries allocated on a TAGE misprediction (1+nnn)

    static const int phistwidth = 27;  // width of the path history used in TAGE
    const int uwidth;  // u counter width on TAGE (2 bits not worth
                           // the effort for a 512 Kbits predictor 0.2 %)
    const int cwidth;  // predictor counter width on the TAGE tagged tables

    // the counter(s) to chose between longest match and alternate prediction on
    // TAGE when weak counters
    unsigned tageConf; // Different confidence levels
    unsigned baseConf;
    unsigned altConf;
    int provVal;
    // For chooser between alt and longest match
    const int alt_width = 5;
    const int size_use_alt;
    int8_t use_alt_on_na[10]; // 10 is not the actual size, but the maximum.

    int TICK;  // for the reset of the u counter

    long long phist;                    // path history

    // History ---------------------------
    const int nHistBits = 2; // Number of history bits per branch
    const bool takenHistory = false; // Use taken history instead of direction history

    // The global history register
    HistoryRegisterFast ghr;
    uint64_t ghr_l;

    // The folded global history registers per table implemented as
    // Circular Shift Registers (CSRs). Each table has three CSRs.
    // One to compute the index and two for the tag. Two are required as any
    // periodic history pattern matching the length of the CSR will XOR to
    // all zero. Therefore the second CSR has a width n-1.
    FoldedHistoryFast* indexFHist[MAXNHIST];
    FoldedHistoryFast* tag1FHist[MAXNHIST];
    FoldedHistoryFast* tag2FHist[MAXNHIST];

    // For the TAGE predictor
    Bentry *btable;             // bimodal TAGE table
    Gentry *gtable[MAXNHIST];  // tagged TAGE tables
    int m[MAXNHIST];
    int TB[MAXNHIST];
    int logg[MAXNHIST];

    int GI[MAXNHIST];  // indexes to the different tables are computed only once

    uint GTAG[MAXNHIST];     // tags for the different tables are computed only once
    bool pred_taken;  // prediction
    bool alttaken;    // alternate  TAGEprediction
    bool tage_pred;   // TAGE prediction
    bool LongestMatchPred;
    int HitBank;  // longest matching bank
    int AltBank;  // alternate matching bank
    Gentry* HitEntry;  // pointer to the HitBank entry
    Gentry* AltEntry;  // pointer to the AltBank entry
    int Seed;     // for the pseudo-random number generator
    bool pred_inter;
    enum {
        BASE    = 0b001,
        ALT     = 0b010,
        LONGEST = 0b100,
        ALT_B   = 0b011,
        ALT_T   = 0b110,
        LAST_TAGE_PROVIDER_TYPE = 0b111
    };

    unsigned provider, tage_provider;

    void reinit();
    int predictorsize();

  public:
    void getPredInfo(unsigned& provider, unsigned& conf) {
        provider = tage_provider;
        conf = (tage_provider == LONGEST) ? tageConf :
               (tage_provider == ALT) ? altConf : baseConf;
    }
  protected:
    // int THRES;

    // Base predictor functions ---------------------------
    // Can be overridden by derived classes
    virtual bool basePredict(const uint64_t pc);
    virtual void baseUpdate(uint64_t pc, bool resolveDir, bool predDir);

    int BI;           // index of the bimodal table
    int8_t BIM;
    bool base_pred;  // prediction of the base table

    // the index functions for the tagged tables uses path history as in the
    // OGEHL predictor
    // F serves to mix path history: not very important impact
    int F(long long A, int size, int bank);

    // gindex computes a full hash of PC, ghist and phist
    int gindex(unsigned int PC, int bank);

    //  tag computation
    uint16_t gtag(unsigned int PC, int bank);

    // Calculate indices and tags for the TAGE predictor
    void calcIndicesAndTags(uint64_t pc);

    // just a simple pseudo random number generator: use available information
    //  to allocate entries  in the loop predictor
    int MYRANDOM();

    // The overall prediction function ---------------------------
    virtual bool predict(uint64_t pc);

    //  TAGE PREDICTION: the actual computation and lookup in the tage tables
    virtual void tagePredict(uint64_t pc);

    // Update the predictor
    virtual void updateTables(uint64_t pc, bool resolveDir, bool predDir);

    // Update of the tagged tables. Will return whether to update the base
    bool tageUpdate(uint64_t pc, bool resolveDir);

    // New table allocations
    virtual int allocate(int idx, uint64_t pc, bool taken);
    virtual int adjustAlloc(bool taken);
    void allocateTables(int nalloc, uint64_t pc, bool taken);
    virtual Gentry& getEntry(int bank);

    int idxChooser();
    virtual unsigned chooseProvider();
    virtual void updateChooser(bool taken);

    // Function to determine when a entry is not useful
    // anymore.
    virtual bool isNotUseful(bool taken);
    virtual bool isUseful(bool taken);


    // History update function
    virtual void updateHistory(const uint64_t pc, const bool taken,
                       const OpType opType, const uint64_t branchTarget);
    virtual void updateGHist(const bool bit);

    virtual bool isAllias(uint64_t pc, int bank);

    virtual void evict(Gentry& entry, int idx) {};

    // Disable bank interleaving.
    const bool disableInterleaving;

    uint64_t branchCount = 0;
    int ticks = 0;

    struct tage_stats {
        int bimodalCorrect = 0;
        int bimodalWrong = 0;
        int longestMatchCorrect = 0;
        int longestMatchWrong = 0;
        int altMatchCorrect = 0;
        int altMatchWrong = 0;
        int total = 0;
        int providerCorrect[MAXNHIST] = {0};
        int providerWrong[MAXNHIST] = {0};
        int allocations[MAXNHIST] = {0};
        int totalAllocations = 0;
        int useful[MAXNHIST] = {0};
        int totalUseful = 0;
        int positiveAlliasing[MAXNHIST] = {0};
        int negativeAlliasing[MAXNHIST] = {0};
        int uResets = 0;
        int totalAlloc = 0;
        int totalAllocInit = 0;
        int totalAllocTries = 0;

        int baseProvTageWouldHaveCorrect = 0;
        int baseProvAltWouldHaveCorrect = 0;
        int baseProvTageSame = 0;
        int baseProvAltSame = 0;
        int altProvTageWouldHaveCorrect = 0;
        int altProvBaseWouldHaveCorrect = 0;
        int altProvTageSame = 0;
        int altProvBaseSame = 0;
        int tageProvBaseWouldHaveCorrect = 0;
        int tageProvAltWouldHaveCorrect = 0;
        int tageProvBaseSame = 0;
        int tageProvAltSame = 0;


        int baseMispred = 0;
        int tageMispred = 0;

        int condBranches = 0;
        int takenBranches = 0;
    } stats;

    virtual void updateStats(bool taken, bool predtaken, uint64_t PC);
};


// Configuration for the TAGE predictor
// Default configuration is for the 64k TAGE predictor
struct TageConfig {

    int nhist = 36;
    int minhist = 6;
    int maxhist = 3000;
    int LogG = 10;
    int LogB = 13;
    int Tbits = 8;
    int nbanklow = 10;
    int nbankhigh = 20;
    int born = 13;
    int assoc_start = 9;
    int assoc_end = 23;

    int uwidth = 1;
    int cwidth = 3;
    int log_size_use_alt = 4;

    bool tage8k = false;
    bool disableInterleaving = false;
    bool overwriteNotUseful = true;
    bool removeAlliasing = false;
    bool tagContext = false;

    void print() const {
        printf("TAGE Config: nhist=%d minhist=%d maxhist=%d LogG=%d LogB=%d Tbits=%d nbanklow=%d nbankhigh=%d born=%d assoc_start=%d assoc_end=%d\n",
            nhist, minhist, maxhist, LogG, LogB, Tbits, nbanklow, nbankhigh, born, assoc_start, assoc_end);
    }
};

inline const TageConfig Tage64kConfig = {};


inline const TageConfig TageLargeConfig = {
    .LogG = 21,
    .LogB = 21,
    .Tbits = 15
};

inline const TageConfig TageLargeBimConfig = {
    .LogB = 21,
};


class Tage64k : public TageBase {
   public:
    Tage64k(void)
        : TageBase(Tage64kConfig) {}
};



};  // namespace LLBP