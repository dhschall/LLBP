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
 */

#pragma once

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <list>
#include <queue>


#include "bpmodels/tage/tage_scl.h"

#include "utils/fileutils.h"
#include "bpmodels/components/counters.h"
#include "bpmodels/components/cache.h"
#include "utils/histogram.h"

namespace LLBP {


struct LLBPConfig;

class LLBP : public TageSCL {

  public:
    LLBP(LLBPConfig config);
    ~LLBP();

    void UpdatePredictor(uint64_t PC, bool resolveDir,
                         bool predDir, uint64_t branchTarget) override;
    void TrackOtherInst(uint64_t PC, OpType opType, bool taken,
                        uint64_t branchTarget) override;
    void PrintStat(double NUMINST) override;
    void tick() override;
    void btbMiss() override;
    void setState(bool warmup) override;

  private:


    // Override some base class functions
    bool predict(uint64_t pc) override;
    void updateTables(uint64_t pc, bool resolveDir, bool predDir) override;
    void updateStats(bool taken, bool predtaken, uint64_t PC) override;
    void updateGHist(const bool bit) override;
    int allocate(int idx, uint64_t pc, bool taken) override;

    bool tageUpdateL1(uint64_t pc, bool resolveDir);
    bool isAllias(uint64_t pc, int bank) override {
        return TageBase::isAllias(pc, bank);
    }

    void resetStats() override;

    void init_predictor();
    void predictorsize();

    int nBrSinceSwitch = 0;
    int mispSinceSwitch = 0;

    typedef uint64_t HistInfo;
    typedef uint64_t Key;
    Key KEY[MAXNHIST];  //

    HistInfo makeHistInfo(int histLength, uint histTag) {
        return ((uint64_t)histTag << 32ULL) | ((uint64_t)histLength);
    }

    struct HistEntry {
      int length;
      uint tag;
      int idx;
      int8_t ctr;
      uint replace;
      bool dir;
      int useful = 0;
      int correct = 0;
      int incorrect = 0;
      Key key = 0;
      int evicted = 0;
      int evicted_ctx = 0;
      int lastEvicted = 0;
      int evictReason = 0;
      uint64_t pc = 0;
      int firstUseDist = -1;
      int maxDist = -1;
    };

    #define MAXHISTSTORE 100


    // Prediction part.
    enum {
      S1 = TageSCL::LAST_SCL_PROVIDER_TYPE + 1,
      S1_BASE,
      S1_TAGE,
      LAST_LCF_PROVIDER_TYPE = S1_TAGE
    };



    // S2 ----------------
    // Second level
    void L2Predict(uint64_t pc);
    void L2Update(uint64_t PC, bool resolveDir, bool predDir);
    bool L2Allocate(int idx, uint64_t pc, bool taken, bool from_l1=false);

    struct Context;


    const int CtrWidth;
    const int ReplCtrWidth;
    const int CtxReplCtrWidth;

    int hitVal[MAXHISTSTORE+1] = {0};
    int HitIdx;  // index of the table that provided the longest matching
                 // prediction
    HistEntry* L2HitEntry;  // pointer to the entry that provided the longest
                          // matching prediction

    // Chooser -----------

    bool isNotUseful(bool taken) override;
    bool isUseful(bool taken) override;
    void updateL2Usefulness(bool taken);

    unsigned chooseProvider() override;


    struct L2PredictorInfo {
        int pVal;
        bool pred;
        unsigned conf;
        int histLength;
        bool promoted;
        bool prefetched;

        bool hit;
        bool isProvider;
        bool shorter;
    } l2;

    bool bim_pred;
    unsigned bimConf;

    // Context

    const int window = 120;
    std::list<uint64_t> bbv[10];

    bool updateRuntimeHash(uint64_t pc, OpType type, bool taken);
    uint64_t calcHash(vector<uint64_t> &vec, int n);
    uint64_t calcHash(std::list<uint64_t> &vec, int n, int start=0, int shift=0);


    uint64_t getCurContext(uint64_t pc);


    ///////////////////////////////////////////////////////
    // LLBP bulk storage


    class PatternSet : public BaseCache<uint64_t, HistEntry>{
    public:

        PatternSet(size_t max_size, size_t assoc) :
            BaseCache<uint64_t, HistEntry>(max_size, assoc)
        {
        }

        HistEntry* insert(const uint64_t &key) {
            return BaseCache<uint64_t, HistEntry>::insert(key);
        }

    };


    // This is the context structure. It contains all the information
    // and statistics for a single context.
    // It also includes the pattern set.
    struct Context {
        bool valid;
        uint64_t key;
        uint64_t pc;
        int correct;
        int incorrect;
        int useful;
        int conflict;
        uint replace;
        int ctr;
        int usefulPtrns;
        int lastPFHit;
        int firstUseDist;
        int maxDist;

        PatternSet patterns;

        Context(uint64_t k, uint64_t p, int n, int assoc)
          : valid(true), key(k), pc(p),
            correct(0), incorrect(0), useful(0), conflict(0),
            replace(0), ctr(0), usefulPtrns(0),
            lastPFHit(0), firstUseDist(-1), maxDist(-1),
            patterns(n, assoc)
        {}

        void sortPatters(const uint64_t key) {
            auto& set = patterns.getSet(key);
            set.sort(
                [](const std::pair<uint64_t, HistEntry>& a, const std::pair<uint64_t, HistEntry>& b)
                {
                    return abs(center(a.second.ctr)) > abs(center(b.second.ctr));
                });
        }

    };



    // The LLBPs large high-capacity structure to store pattern sets.
    // Its implemented as a set associative cache.
    // The Context directory (CD) can be thought of as the tag array while the
    // LLBPStorage is the data array. In the implementation there is only
    // one data structure.
    class LLBPStorage : public BaseCache<uint64_t, Context>{
        typedef typename std::pair<uint64_t, Context> key_value_pair_t;
	    typedef typename std::list<key_value_pair_t>::iterator list_iterator_t;
        const int n_patterns;
        const int _ptrn_assoc;

    public:

        LLBPStorage(int n_ctx, int n_patterns, int ctx_assoc, int ptrn_assoc)
          : BaseCache<uint64_t, Context>(n_ctx, ctx_assoc),
            n_patterns(n_patterns), _ptrn_assoc(ptrn_assoc)
        {
        }

        // Will create a new context but does not install it
        Context* createNew(uint64_t key, uint64_t pc) {
            return new Context(key, pc, n_patterns, _ptrn_assoc);
        }

        Context* insert(uint64_t key, uint64_t pc) {

		    auto c = this->get(key);
            if (c != nullptr) {
                return c;
            }

            auto& set = this->getResizedSet(key);

            set.push_front(
                key_value_pair_t(key, Context(key, pc, n_patterns, _ptrn_assoc)));
            _index[key] = set.begin();
            return &set.front().second;
        }

        void sortContexts(uint64_t key) {
            auto& set = this->getSet(key);
            set.sort(
                [](const key_value_pair_t& a, const key_value_pair_t& b)
                {
                    return a.second.replace > b.second.replace;
                });
        }
    };

    LLBPStorage llbpStorage;

    Context* HitContext;

    const int numContexts;
    const int numPatterns;


    int lastSwitch = 0;

    // The variables to define the context hash function
    const int hashT, hashN, hashSk, hashSh;

    struct CTX_t {
        uint64_t cur = 0;
        uint64_t prev = 0;
        uint64_t lastSwitch = 0;
        uint64_t slackSamples = 0;
        uint64_t totalSlack = 0;
        uint64_t lastUseful = 0;
        uint64_t prevUseful = 0;
        uint64_t usefulSamples = 0;
        uint64_t totalUseful = 0;
    };

    CTX_t contexts[80] = {0};

    void updateContext(CTX_t &ctx, uint64_t cur);

    Context* allocateNewContext(uint64_t pc, uint64_t key);

    // PCHistoryRegister* pcHistory;

    Key lastHitCtx = 0;
    int lastCtxSwitch = 0;

    std::unordered_map<uint64_t,std::unordered_map<uint64_t,HistEntry>> allPatterns;
    void evictPattern(uint64_t ctx_key, HistEntry* ptrn, bool ctx_evict=false);

    // We do not use all tables (history lengths) in LLBP
    struct PBEntry {
        Key key;
        bool dirty;
        bool newlyAllocated;
        bool used;
        bool useful;
        int origin;
        bool valid;
        int prefetchtime;
        bool locked;
        PBEntry(Key c)
          : key(c), dirty(false),
            newlyAllocated(false),
            used(false), useful(false),
            origin(0),
            valid(false),
            prefetchtime(0),
            locked(false)
        {}
        PBEntry() : PBEntry(0) {}
    };
    PBEntry* L2PredHit;
    class PatternBuffer : public BaseCache<uint64_t, PBEntry> {
      public:
        PatternBuffer(int n, int assoc)
          : BaseCache<uint64_t, PBEntry>(n, assoc)
        {
        }

        PBEntry* insert(PBEntry &entry) {;
            auto v = get(entry.key);
            if (v != nullptr) {
                return v;
            }

            // Get the set with a free item
            auto& set = getResizedSet(entry.key);

            set.push_front(key_value_pair_t(entry.key, entry));
            _index[entry.key] = set.begin();
            return &set.front().second;
        }
    };

    PatternBuffer patternBuffer;

    void prefetch();
    void tickPrefetchQueue();
    void squashPrefetchQueue(bool btbMiss=false);
    void installInPB(PBEntry &entry, bool bypass=false);
    std::list<PBEntry> prefetchQueue;


    const bool simulateTiming;
    bool warmup = false;
    const int accessDelay;

    int lastMispredict = 0;

    std::unordered_map<int,int> fltTables;



    const int TTWidth;
    const int CTWidth;
    FoldedHistoryFast* fghrT1[MAXNHIST];
    FoldedHistoryFast* fghrT2[MAXNHIST];


    // Some histograms
    Histogram<size_t> primMispLength;
    Histogram<size_t> llbpMispLength;
    Histogram<size_t> primProvLength;
    Histogram<size_t> llbpProvLength;
    Histogram<size_t> numHistPerContext;
    Histogram<size_t> numUsefulHistPerContext;


    // Some stats
    struct l2_stats
    {

        int tageProv = 0;
        int tageCorrect = 0;
        int tageWrong = 0;
        int sclProv = 0;
        int sclCorrect = 0;
        int sclWrong = 0;
        int baseProv = 0;
        int baseCorrect = 0;
        int baseWrong = 0;



        // int override = 0;
        // int goodOverride = 0;
        // int badOverride = 0;
        // int badOverrideLonger = 0;
        // int badOverrideSame = 0;
        // int l2Hit = 0;
        // int overrideWhouldHaveCorrect = 0;
        // int overrideWhouldHaveWrong = 0;
        // int bim1WhouldHaveCorrect = 0;
        // int bim1WhouldHaveWrong = 0;
        // int bim2WhouldHaveCorrect = 0;
        // int bim2WhouldHaveWrong = 0;
        // int bimContMatch = 0;
        // int overrideSame = 0;
        // int overrideNoHit = 0;
        // int clearAlloc = 0;
        // int clearUbit = 0;
        // int collision = 0;
        // int baseMispred = 0;
        // int tageMispred = 0;
        // int s1Mispred = 0;
        // int sclMispred = 0;

        int l2Prov = 0;
        int l2Correct = 0;
        int l2Wrong = 0;
        int l2CtxHit = 0;
        int l2PtrnHit = 0;
        // int l2Alloc = 0;
        // int l2exists = 0;
        // int l1useful = 0;
        // int l1useless = 0;
        // int l1bimuseful = 0;
        // int l1bimuseless = 0;
        // int l2useful = 0;
        // int l2useless = 0;
        // int l1l2correct = 0;
        // int l1l2hitOnSame = 0;
        // int l2longer = 0;
        // int l2uselessSignChange = 0;
        // int l2uselessAltSame = 0;
        // int l2uselessBimSame = 0;
        // int l2notBecauseNotUseful = 0;
        int l2notBecauseShorter = 0;
        // int l2usefulButNotPromoted = 0;
        int l2notBecauseNotPrefetched = 0;
        // int noAllocBecauseL2CorrNoPromote = 0;

        // int l2CorrectL1Wrong = 0;
        // int l1NoHit = 0;
        // int l1EvictDueToSignChange = 0;
        // int l1EvictDueToAltSame = 0;
        // int l1EvictDueToL2Same = 0;
        // int l1EvictBeforeUseful = 0;
        // int l1EvictUReset = 0;
        // int UseEvictDueToSignChange = 0;
        // int UseEvictDueToAltSame = 0;
        // int UseEvictDueToL2Same = 0;
        // int UseEvictBeforeUseful = 0;
        // int UseEvictUReset = 0;
        // int l1HitSameButWrong = 0;
        // int l1HitSameCtrDiff = 0;
        // int l1HitSameUseAlt = 0;
        // int l1promote = 0;
        // int l1promoteSuccess = 0;
        // int l1promoteFail = 0;
        // int l1promoteExits = 0;

        // int ctxSwitches = 0;
        // int ctxSame = 0;
        // int ctxCached = 0;
        // int ctxSwitchesUseful = 0;
        // int ctxSwitches1Useful = 0;
        // uint64_t numPatternsFetched = 0;
        // uint64_t numUPatternsFetched = 0;
        // uint64_t numPatternsFetchedUseful = 0;
        // uint64_t numUPatternsFetchedUseful = 0;

        // int belowInsert = 0;
        // int aboveInsert = 0;

        // int numRelearn = 0;
        // int wrongInit = 0;
        // int evictSignChangeB1 = 0;
        // int evictAltSameB1 = 0;

        // int l2ctxExists = 0;
        // int l2ctxNoExists = 0;
        // int patternBufferd = 0;
        // int l2predInL2 = 0;
        int l2cacheDirtyEvict = 0;
        int l2cacheCleanEvict = 0;
        // int l2cacheNewAlloc = 0;
        // int l2cacheUsefulEvict = 0;
        // int l2cacheUsed = 0;

        // int l2evicts = 0;
        // int l2evictsUseful = 0;
        // int l2evictNumPatterns = 0;
        // int l2evictFromCache = 0;
        // int l2evictFromQueue = 0;

        int l2PFHitInQueue = 0;
        int l2PFHitInCache = 0;
        int l2PFHitInCI = 0;
        int pfDroppedLocked = 0;
        int pfDroppedMispredict = 0;
        int pfDroppedBTBMiss = 0;

        int primCorrect = 0;
        int primWrong = 0;
        int l2OverrideGood = 0;
        int l2OverrideBad = 0;
        int l2OverrideSameCorr = 0;
        int l2OverrideSameWrong = 0;
        int l2NoOverride = 0;
        int ovrPosAlias = 0;
        int ovrNegAlias = 0;

        // int fullCtxAlloc = 0;
        // int partCtxAlloc = 0;
        // int fullCtxUseful = 0;
        // int partCtxUseful = 0;

    } l2stats;

};


struct LLBPConfig {
  TSCLConfig tsclConfig;


#define LLBP_CONSTRAINED
#ifdef LLBP_CONSTRAINED
  int numPatterns = 16;
#else
  int numPatterns = 1000000;
#endif


#ifdef LLBP_CONSTRAINED
  int numContexts = 1024*14;
#else
  int numContexts = 1000000;
#endif


#ifdef LLBP_CONSTRAINED
    int ctxAssoc = 7;
    int ptrnAssoc = 4;
#else
    int ctxAssoc = numContexts;
    int ptrnAssoc = numPatterns;
#endif

    int CtrWidth = 3;
    int ReplCtrWidth = 16;
    int CtxReplCtrWidth = 2;

    int l2cacheSize = 64;

#ifdef LLBP_CONSTRAINED
    int l2cacheAssoc = 4;
#else
    int l2cacheAssoc = l2cacheSize;
#endif

#ifdef LLBP_CONSTRAINED
    int TTWidth = 13;
    int CTWidth = 14;
#else
    int TTWidth = 20;
    int CTWidth = 31;
#endif

    bool pcInCtx = false;
    bool simulateTiming = false;
    int accessDelay = 6;


  void print() const {
    printf("LLBP Config: NumPatterns=%i, NumContexts=%i, ctxAssoc=%i, ptrnAssoc=%i, CtrWidth=%i, ReplCtrWidth=%i, CtxReplCtrWidth=%i, l2cacheSize=%i, TTWidth=%i, CTWidth=%i, simMispFlush=%i, accessDelay=%i\n ",
           numPatterns, numContexts, ctxAssoc, ptrnAssoc, CtrWidth, ReplCtrWidth, CtxReplCtrWidth, l2cacheSize, TTWidth, CTWidth, simulateTiming, accessDelay);
  }
};




/////////////////////////
// LLBP Predictor

// The LLBP predictor without simulating the prefetch latency.
class LLBPTageSCL64k : public LLBP {
   public:
    LLBPTageSCL64k(void)
        : LLBP(LLBPConfig
            {
                .tsclConfig = TSCLConfig
                {
                    .tageConfig = Tage64kConfig,
                    .useSC = true,
                    .useLoop = true
                },
                .simulateTiming = false,
            })
    {}
};


// The LLBP predictor with simulating the prefetch latency.
class LLBPTageSCL64kTiming : public LLBP {
   public:
    LLBPTageSCL64kTiming(void)
        : LLBP(LLBPConfig
            {
                .tsclConfig = TSCLConfig
                {
                    .tageConfig = Tage64kConfig,
                    .useSC = true,
                    .useLoop = true
                },
                .simulateTiming = true,
            })
    {}
};







};  // namespace LLBP
