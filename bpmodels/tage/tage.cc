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

#include "tage.h"

#include "utils/common.h"
#include "bpmodels/components/counters.h"

using namespace std;

#define BORNTICK 1024
#define PRINTDEBUG 0

#define COND (stats.total >= 00000)

#define USEPATH

namespace LLBP {


TageBase::TageBase(TageConfig cfg)
    : BasePredictor(),
        nhist(cfg.nhist),
        nbanklow(cfg.nbanklow),
        nbankhigh(cfg.nbankhigh),
        born(cfg.born),
        assoc_start(cfg.assoc_start),
        assoc_end(cfg.assoc_end),
        minhist(cfg.minhist),
        maxhist(cfg.maxhist),
        LogB(cfg.LogB),
        LogG(cfg.LogG),
        Tbits(cfg.Tbits),
        uwidth(cfg.uwidth),
        cwidth(cfg.cwidth),
        size_use_alt(1 << (cfg.log_size_use_alt)),
        ghr(histbufferlen),
        disableInterleaving(cfg.disableInterleaving)
{
    assert(minhist <= maxhist);
    assert(LogG > 0);
    assert(LogB > 0);
    assert(Tbits > 0);
    assert(nhist <= MAXNHIST);

    cfg.print();

    // initialize the predictor
    reinit();
    predictorsize();
}

TageBase::~TageBase() {
    delete[] btable;

    delete[] gtable[1];
    delete[] gtable[born];

    for (int i = 1; i <= nhist; i++) {
        delete indexFHist[i];
        delete[] tag1FHist[i];
        delete[] tag2FHist[i];
    }
}

void TageBase::reinit() {
    m[1] = minhist;
    m[nhist / 2] = maxhist;
    for (int i = 2; i <= nhist / 2; i++) {
        m[i] = (int)(((double)minhist *
                      pow((double)(maxhist) / (double)minhist,
                          (double)(i - 1) / (double)(((nhist / 2) - 1)))) +
                     0.5);
        //      printf("(%d %d)", m[i],i);
    }
    for (int i = 1; i <= nhist; i++) {
        NOSKIP[i] = ((i - 1) & 1) || ((i >= assoc_start) & (i < assoc_end));
    }

    if (nhist > 30) {
        NOSKIP[4] = 0;
        NOSKIP[nhist - 2] = 0;
        NOSKIP[8] = 0;
        NOSKIP[nhist - 6] = 0;
        // just eliminate some extra tables (very very marginal)
    }

    for (int i = nhist; i > 1; i--) {
        m[i] = m[(i + 1) / 2];
    }
    for (int i = 1; i <= nhist; i++) {
        TB[i] = Tbits + 4 * (i >= born);
        logg[i] = LogG;
    }

    gtable[1] = new Gentry[nbanklow * (1 << LogG)];
    SizeTable[1] = nbanklow * (1 << LogG);

    gtable[born] = new Gentry[nbankhigh * (1 << LogG)];
    SizeTable[born] = nbankhigh * (1 << LogG);

    for (int i = born + 1; i <= nhist; i++) gtable[i] = gtable[born];
    for (int i = 2; i <= born - 1; i++) gtable[i] = gtable[1];


    btable = new Bentry[1 << LogB];

    for (int i = 1; i <= nhist; i++) {
        indexFHist[i] = new FoldedHistoryFast(ghr, m[i], logg[i]);
        tag1FHist[i] = new FoldedHistoryFast(ghr, m[i], TB[i]);
        tag2FHist[i] = new FoldedHistoryFast(ghr, m[i], TB[i] - 1);
    }





    Seed = 0;

    TICK = 0;
    phist = 0;
    Seed = 0;

    // for (int i = 0; i < histbufferlen; i++) ghist[0] = 0;
    // ptghist = 0;
    // updatethreshold = 35 << 3;

    for (int i = 0; i < (1 << LogB); i++) {
        btable[i].pred = 0;
        btable[i].hyst = 1;
    }
    for (int i = 0; i < size_use_alt; i++) {
        use_alt_on_na[i] = 0;
    }

    // ptghist = 0;
    phist = 0;
}

int TageBase::predictorsize() {
    int STORAGESIZE = 0;

    STORAGESIZE +=
        nbankhigh * (1 << (logg[born])) * (cwidth + uwidth + TB[born]);
    STORAGESIZE += nbanklow * (1 << (logg[1])) * (cwidth + uwidth + TB[1]);

    STORAGESIZE += (size_use_alt)*alt_width;
    STORAGESIZE += (1 << LogB) + (1 << (LogB - hystshift));
    STORAGESIZE += m[nhist];
    STORAGESIZE += phistwidth;
    STORAGESIZE += 10;  // the TICK counter

    printf("LogG:%i, TBITS:%i, UWIDTH:%i, CWIDTH:%i, ALTWIDTH:%i, LogB:%i, Hyst:%i\n",
            LogG, Tbits, uwidth, cwidth, alt_width, LogB, hystshift);


    printf(" (TAGE %d) ", STORAGESIZE);

    printf(" (TOTAL %d bits %d Kbits) ", STORAGESIZE,
            STORAGESIZE / 1024);

    // for printing predictor characteristics
    int NBENTRY = 0;

    STORAGESIZE = 0;
    for (int i = 1; i <= nhist; i++) {
        if (NOSKIP[i]) {
            printf("%dx%d ", TB[i], (1 << logg[i]));

            STORAGESIZE += (1 << logg[i]) * (5 + TB[i]);
            NBENTRY += (1 << logg[i]);
        }
    }

    printf("\n");

    for (int i = 1; i <= nhist; i++) {
        if (NOSKIP[i]) printf("%d ", m[i]);
    }
    printf("\n");

    printf("TAGE: N:%d -> SIZE:%d,%iK\n", NBENTRY, STORAGESIZE,
            STORAGESIZE / 1024);

    int BIMSIZE = (1 << LogB) + (1 << (LogB - hystshift));
    printf("BASE: Dir:%i, Hyst:%i -> SIZE: %d, %dK\n", (1 << LogB),
            (1 << (LogB - hystshift)), BIMSIZE, BIMSIZE / 1024);
    STORAGESIZE += BIMSIZE;

    printf("nhist= %d; MInhist= %d; MAXHIST= %d; STORAGESIZE= %d, %dKB; "
            "NBENTRY= %d\n",
            nhist, minhist, maxhist, STORAGESIZE, STORAGESIZE / 1024, NBENTRY);

    return (STORAGESIZE);
}

// Base Predictions ----------------------------------------
bool TageBase::basePredict(const uint64_t pc) {
    //
    tage_provider = BASE;
    BI = (pc ^ (pc >> 2)) & ((1 << LogB) - 1);

    provVal = BIM = (btable[BI].pred << 1) + (btable[BI >> hystshift].hyst);
    baseConf = ((BIM == 0) || (BIM == 3)) ? HighConf : LowConf;
    base_pred = BIM > 1;

    return base_pred;
}

void TageBase::baseUpdate(uint64_t pc, bool resolveDir, bool predDir) {
    int inter = BIM;
    if (resolveDir) {
        if (inter < 3) inter += 1;
    } else if (inter > 0)
        inter--;
    btable[BI].pred = inter >> 1;
    btable[BI >> hystshift].hyst = (inter & 1);

    // ctrupdate(btable[BI].ctr, resolveDir, 2);
    btable[BI].pc = pc;
}


int TageBase::F(long long A, int size, int bank) {
    int A1, A2;
    A = A & ((1 << size) - 1);
    A1 = (A & ((1 << logg[bank]) - 1));
    A2 = (A >> logg[bank]);

    if (bank < logg[bank])
        A2 = ((A2 << bank) & ((1 << logg[bank]) - 1)) +
             (A2 >> (logg[bank] - bank));
    A = A1 ^ A2;
    if (bank < logg[bank])
        A = ((A << bank) & ((1 << logg[bank]) - 1)) +
            (A >> (logg[bank] - bank));
    return (A);
}

// gindex computes a full hash of PC, ghist and phist
int TageBase::gindex(unsigned int PC, int bank) {
    int index;
    index = PC ^ (PC >> (abs(logg[bank] - bank) + 1));

    index ^= indexFHist[bank]->value;


#ifdef USEPATH
    int M = (m[bank] > phistwidth) ? phistwidth : m[bank];
    index ^= F(phist, M, bank);
#endif

    return (index & ((1 << (logg[bank])) - 1));
}

//  tag computation
uint16_t TageBase::gtag(unsigned int PC, int bank) {

    int tag = 0;
    tag = PC;
    tag ^= tag1FHist[bank]->value ^ (tag2FHist[bank]->value << 1);

    return (tag & ((1 << (TB[bank])) - 1));
}

// just a simple pseudo random number generator: use available information
//  to allocate entries  in the loop predictor
int TageBase::MYRANDOM() {
    Seed++;
    Seed ^= phist;
    Seed = (Seed >> 21) + (Seed << 11);
    Seed = (Seed >> 10) + (Seed << 22);
    return (Seed);
};



// Prediction ----------------------------------------------
bool TageBase::GetPrediction(uint64_t PC) {

    DPRINTIF(COND,"---- %i PC:%lx -------\n", stats.total, PC);
    // computes the TAGE table addresses and the partial tags
    pred_taken = predict(PC);
    return pred_taken;
}

bool TageBase::predict(uint64_t PC) {

    // 1. The base prediction
    basePredict(PC);
    tage_provider = BASE;

    // 2. The TAGE prediction
    tagePredict(PC);
    provider = tage_provider;
    return tage_pred;
}


void TageBase::calcIndicesAndTags(uint64_t PC) {

    // 1. Compute indices and tags
    for (int i = 1; i <= nhist; i += 2) {
        GI[i] = gindex(PC, i);
        GTAG[i] = gtag(PC, i);
        GTAG[i + 1] = GTAG[i];
        GI[i + 1] = GI[i] ^ (GTAG[i] & ((1 << LogG) - 1));
    }


    int T = (PC ^ (phist & ((1 << m[born]) - 1))) % nbankhigh;
    T = disableInterleaving ? 1 : T;

    for (int i = born; i <= nhist; i++)
        if (NOSKIP[i]) {
            GI[i] += (T << LogG);
            T++;
            T = T % nbankhigh;
        }

    T = (PC ^ (phist & ((1 << m[1]) - 1))) % nbanklow;
    T = disableInterleaving ? 1 : T;

    for (int i = 1; i <= born - 1; i++)
        if (NOSKIP[i]) {
            GI[i] += (T << LogG);
            T++;
            T = T % nbanklow;
        }

}

//  TAGE PREDICTION: same code at fetch or retire time but the index and
//  tags must recomputed
void TageBase::tagePredict(uint64_t PC) {

    // 1. Compute indices and tags
    calcIndicesAndTags(PC);

    HitBank = AltBank = 0;
    HitEntry = AltEntry = nullptr;
    #define TAGBRANCH (PC == 94321081556116)

    // 2. Perform the table lookup
    // Look for the bank with longest matching history
    for (int i = nhist; i > 0; i--) {
        if (NOSKIP[i]) {
            if (gtable[i][GI[i]].tag == GTAG[i]) {

                HitBank = i;
                HitEntry = &gtable[i][GI[i]];
                break;
            }
        }
    }

    // Look for the alternate bank
    for (int i = HitBank - 1; i > 0; i--) {
        if (NOSKIP[i]) {
            if (gtable[i][GI[i]].tag == GTAG[i]) {

                AltBank = i;
                AltEntry = &gtable[i][GI[i]];
                break;
            }
        }
    }

    // If there was no hit in the tables use the base prediction
    // Initialize the prediction to the base prediction
    tageConf = altConf = baseConf;
    LongestMatchPred = tage_pred = alttaken = base_pred;

    // 3. Read the predictions and choose between
    // longest matching and alternate matching
    if (HitBank > 0) {

        // Read the longest match prediction and its confidence.
        LongestMatchPred = (HitEntry->ctr >= 0);
        tageConf = compConf(HitEntry->ctr, cwidth);

        if (AltBank > 0) {
            // For a second hit read also the alternate match prediction.
            alttaken = (AltEntry->ctr >= 0);
            altConf = compConf(AltEntry->ctr, cwidth);
        }
        DPRINTIF(COND,"Hit:%i,GI:%i,GT:%i,c:%i Alt:%i\n", HitBank, GI[HitBank], GTAG[HitBank], HitEntry->ctr, AltBank);


        // Manage the selection between longest matching and alternate
        // matching is done by considering the confidence of longest
        // and alternate matching.
        tage_provider = chooseProvider();
        switch (tage_provider) {
            case LONGEST:
                provVal = HitEntry->ctr;
                tage_pred = LongestMatchPred;
                break;
            case ALT:
                provVal = AltEntry->ctr;
                tage_pred = alttaken;
                break;
            case BASE:
                provVal = BIM;
                tage_pred = base_pred;
                break;
        }
    }
}


void TageBase::updateHistory(const uint64_t pc, const bool taken,
                                const OpType opType, const uint64_t target) {

    bool indirect = (opType == OPTYPE_CALL_INDIRECT_COND) ||
                    (opType == OPTYPE_CALL_INDIRECT_UNCOND) ||
                    (opType == OPTYPE_JMP_INDIRECT_COND) ||
                    (opType == OPTYPE_JMP_INDIRECT_UNCOND) ||
                    (opType == OPTYPE_RET_COND) ||
                    (opType == OPTYPE_RET_UNCOND);

    int tbits = indirect ? nHistBits+1 : nHistBits;

    int historyBits = taken ? 1 : 0;
    if (nHistBits > 1) {
        historyBits ^= pc ^ (pc >> 2);
    }
    ghr_l = (ghr_l << 1) | (taken & 1);
    int PATH = pc ^ (pc >> 2) ^ (pc >> 4);

    for (int t = 0; t < tbits; t++) {
        // update  history

        bool brDir = (historyBits & 1);
        updateGHist(brDir);
        historyBits >>= 1;
#ifdef USEPATH
        int PATHBIT = (PATH & 127);
        PATH >>= 1;
        phist = (phist << 1) ^ PATHBIT;
        phist &= (1<<phistwidth)-1;
#endif
    }
}


void TageBase::updateGHist(const bool bit) {
    ghr.push(bit);
    for (int i = 1; i <= nhist; ++i) {
        indexFHist[i]->update();
        tag1FHist[i]->update();
        tag2FHist[i]->update();
    }
}


int TageBase::idxChooser() {
    bool add1 = (altConf != LowConf);
    return ((((HitBank - 1) / 8) << 1) + add1) % (size_use_alt - 1);
}

unsigned TageBase::chooseProvider() {
    // Manage the selection between longest matching and alternate
    // matching.
    // We take two sources of information to decide. First wheather
    // the hit entry is recognized as a newly allocated entry (confidence
    // low) and USE_ALT_ON_NA is positive  use the alternate prediction

    // If the longest is somehow certain use its prediction.
    if (tageConf != LowConf) {
        return LONGEST;
    }

    // Use on low confidence if the USE_ALT_ON_NA is negative
    if (use_alt_on_na[idxChooser()] < 0) {
        return LONGEST;
    }

    return (AltBank > 0) ? ALT : BASE;
}

void TageBase::updateChooser(bool taken) {

    if (tageConf != LowConf)
        return;

    if (LongestMatchPred != alttaken) {
        ctrupdate(use_alt_on_na[idxChooser()], (alttaken == taken), alt_width);
    }
}


int
TageBase::adjustAlloc(bool resolveDir) {

    // TAGE UPDATE
    bool ALLOC = ((tage_pred != resolveDir) & (HitBank < nhist));

    // do not allocate too often if the overall prediction is correct

    if (HitBank > 0) {
        // Manage the selection between longest matching and alternate
        // matching for "pseudo"-newly allocated longest matching entry this
        // is extremely important for TAGE only, not that important when the
        // overall predictor is implemented
        // An entry is considered as newly allocated if its prediction
        // counter is weak
        if (tageConf == LowConf) {

            // if the longest match was delivering the correct prediction,
            // no need to allocate a new entry even if the overall
            // prediction was false
            if (LongestMatchPred == resolveDir) ALLOC = false;
        }
        // Update the chooser policy between longest matching and
        // alternate matching.
        updateChooser(resolveDir);
    }

    if (tage_pred == resolveDir)
        if ((MYRANDOM() & 31) != 0) ALLOC = false;

    // Fixed number of allocations for a missprediction.
    if (ALLOC) {
        return 1 + nnn;
    }
    return 0;
}


void
TageBase::allocateTables(int nalloc, uint64_t PC, bool resolveDir) {
    if (nalloc <= 0) return;

    // T is the number of entries additionally allocated to at
    // least one entry per missprediction.
    int T = nalloc - 1;

    int A = 1;
    if ((MYRANDOM() & 127) < 32) A = 2;
    int Penalty = 0;
    int NA = 0;
    // int
    int DEP = HitBank + 1;
    if (assoc_start < assoc_end)
        DEP = ((((HitBank - 1 + 2 * A) & 0xffe)) ^ (MYRANDOM() & 1));

    // just a complex formula to chose between X and X+1, when X is odd:
    // sorry

    for (int I = DEP; I < nhist; I += 2) {
        int i = I + 1;
        bool Done = false;
        if (NOSKIP[i]) {
            auto n = allocate(i, PC, resolveDir);
            if (n > 0) {
                NA+=1;

                if ((T <= 0) || n > 1) {
                    break;
                }
                I += 2;
                Done = true;
                T -= 1;
            } else if (n < 0) {
                Penalty++;
            }
        }

        if (!Done) {
            i = (I ^ 1) + 1;
            if (NOSKIP[i]) {
                auto n = allocate(i, PC, resolveDir);
                if (n > 0) {
                    NA+=1;

                    if ((T <= 0) || n > 1) {
                        break;
                    }
                    I += 2;
                    Done = true;
                    T -= 1;
                } else if (n < 0) {
                    Penalty++;
                }
            }
        }
    }
    stats.totalAllocTries += Penalty + NA;
    stats.totalAllocInit++;
    TICK += (Penalty - 2 * NA);

    // just the best formula for the Championship:
    // In practice when one out of two entries are useful
    if (TICK < 0) TICK = 0;
    if (TICK >= BORNTICK) {
        for (int i = 1; i <= born; i += born - 1) // born=11 => T:1,11
            for (int j = 0; j < SizeTable[i]; j++) {
                gtable[i][j].u >>= 1;
            }
        TICK = 0;
        stats.uResets++;
    }
}

Gentry& TageBase::getEntry(int idx) {
    return gtable[idx][GI[idx]];
}

int TageBase::allocate(int idx, uint64_t pc, bool taken) {
    auto& entry = getEntry(idx);

    if (entry.u != 0) {
        return -1;
    }


#define OPTREMP
// the replacement is optimized with a single u bit: 0.2 %
#ifdef OPTREMP
    if (abs(2 * entry.ctr + 1) > 3) {
        if (entry.ctr > 0)
            entry.ctr--;
        else
            entry.ctr++;
        return 0;
    }
#endif

    evict(entry, idx);

    DPRINTIF(COND,"Alloc:%i,GI:%i,GT:%i\n", idx, GI[idx], GTAG[idx]);

    entry.tag = GTAG[idx];
    entry.pc = pc;
    entry.hlen = idx;
    entry.idx = GI[idx];

    entry.ctr = (taken) ? 0 : -1;

    // entry.u = 0;
    stats.allocations[idx]++;
    stats.totalAllocations++;

    if (entry.correct < 0) entry.correct = 0;

    return 1;
}


// PREDICTOR UPDATE

bool TageBase::tageUpdate(uint64_t pc, bool resolveDir) {

    bool update_base = false;

    // update predictions
    if (HitBank > 0) {
        if (tageConf == LowConf) {
            if (LongestMatchPred != resolveDir) {
                // acts as a protection
                if (AltBank > 0) {
                    ctrupdate(AltEntry->ctr, resolveDir, cwidth);
                } else {
                    update_base = true;
                }
            }
        }

        // Do the actual counter update
        ctrupdate(HitEntry->ctr, resolveDir, cwidth);
        // sign changes: no way it can have been useful
        if (HitEntry->ctr == (resolveDir ? 0 : -1)) {
            HitEntry->u = 0;
        }

        // If both the longest and alternate predictions where correct
        // we can possible free the longest entry to use it for other
        // predictions.
        // We clear this entry by clearing the useful bit.
        if (isNotUseful(resolveDir)) {
            if (HitEntry->u > 0) {
                HitEntry->u--;
            }
        }
        // If the longest hit was correct but the alternative prediction
        // was not promote this entry to be useful.
        if (isUseful(resolveDir)) {
            if (HitEntry->u < (1 << uwidth) - 1)
                HitEntry->u++;
            HitEntry->useful++;
        }
        DPRINTIF(COND,"TageUpdate: idx:%d, ctr:%i,u:%i,T:%i\n",
                HitBank, HitEntry->ctr, HitEntry->u, HitEntry->tag);

    } else {
        update_base = true;
    }

    // END TAGE UPDATE
    return update_base;
}

bool TageBase::isNotUseful(bool taken) {

    // If both the longest and alternate predictions where correct
    // we can possible free the longest entry to use it for other
    // predictions.
    if ((alttaken == taken) && (LongestMatchPred == taken)) {
        // We only clear if the alternate prediction has a
        // high confidence.
        if (altConf == HighConf) {
            if (AltBank > 0) {
                return true;
            }
        }
    }
    return false;
}

bool TageBase::isUseful(bool taken) {

    // If the longest prediction is correct but the alternate
    // prediction is wrong the longest is useful.
    if ((alttaken != taken) && (LongestMatchPred == taken)) {
        return true;
    }
    return false;
}


void TageBase::updateTables(uint64_t pc, bool resolveDir, bool predDir) {

    // 1. Allocate tables if necessary
    int nalloc = adjustAlloc(resolveDir);
    allocateTables(nalloc, pc, resolveDir);

    // 2. the TAGE tables
    bool update_base = tageUpdate(pc, resolveDir);

    // If the prediction was from the base predictor, update it.
    if (update_base) {
        baseUpdate(pc, resolveDir, predDir);
    }
}

void TageBase::FirstTimeUpdate(uint64_t PC, bool taken,
                                uint64_t branchTarget) {

    branchCount++;
    if (taken) {
        stats.takenBranches++;
    }
    // computes the TAGE table addresses and the partial tags
    basePredict(PC);
    baseUpdate(PC, taken, false);
    updateHistory(PC, taken, OPTYPE_JMP_DIRECT_COND, branchTarget);
}


// Predictor update ----------------------------------------
void TageBase::UpdatePredictor(uint64_t PC, bool resolveDir,
                                  bool predDir, uint64_t branchTarget) {
    branchCount++;
    stats.condBranches++;
    if (predDir) {
        stats.takenBranches++;
    }
    updateTables(PC, resolveDir, predDir);

    updateHistory(PC, resolveDir, OPTYPE_JMP_DIRECT_COND, branchTarget);
    updateStats(resolveDir, predDir, PC);
}

void TageBase::TrackOtherInst(uint64_t PC, OpType opType, bool taken,
                                 uint64_t branchTarget) {
    branchCount++;
    if (taken) {
        stats.takenBranches++;
    }
    updateHistory(PC, taken, opType, branchTarget);
}





bool TageBase::isAllias(uint64_t pc, int bank) {
    if (bank == 0) {
        return (btable[BI].pc != pc) && (btable[BI].pc != 0);
    };
    return (gtable[bank][GI[bank]].pc != pc) && (gtable[bank][GI[bank]].pc != 0);
}

void TageBase::updateStats(bool taken, bool predtaken, uint64_t PC) {
    stats.total++;
    auto correct = taken == tage_pred;
    if (correct) {
        auto cb = 0;
        switch (tage_provider) {
            case LONGEST:
                stats.longestMatchCorrect++;
                cb = HitBank;
                break;
            case ALT:
                stats.altMatchCorrect++;
                cb = AltBank;
                break;
            case BASE:
                stats.bimodalCorrect++;
                cb = 0;
                break;
        }
        stats.providerCorrect[cb]++;
        if (isAllias(PC, cb)) {
            stats.positiveAlliasing[cb]++;
        }
    } else {
        auto cb = 0;
        switch (tage_provider) {
            case LONGEST:
                stats.longestMatchWrong++;
                cb = HitBank;
                break;
            case ALT:
                stats.altMatchWrong++;
                cb = AltBank;
                break;
            case BASE:
                stats.bimodalWrong++;
                cb = 0;
                break;
        }
        stats.providerWrong[cb]++;
        if (isAllias(PC, cb)) {
            stats.negativeAlliasing[cb]++;
        }
    }

    // If it was new allocation and it was delivering the correct
    // prediction update stats
    if (correct) {
        switch (tage_provider) {
            case LONGEST:
                HitEntry->correct++;
                break;

            case ALT:
                AltEntry->correct++;
                break;

            case BASE:
                break;
        }
    } else {
        switch (tage_provider) {
            case LONGEST:
                HitEntry->incorrect++;
                break;

            case ALT:
                AltEntry->incorrect++;
                break;

            case BASE:
                break;
        }
    }

    auto tage_correct = taken == tage_pred;
    auto base_correct = taken == base_pred;
    auto alt_correct = taken == alttaken;
    auto long_correct = taken == LongestMatchPred;

    if (HitBank > 0) {
        if (AltBank > 0) {
            if (tage_provider == ALT) {

                if (alttaken == LongestMatchPred) {
                    stats.altProvTageSame++;
                }
                if (alttaken == base_pred) {
                    stats.altProvBaseSame++;
                }


                if (!alt_correct && long_correct) {
                    stats.altProvTageWouldHaveCorrect++;
                }
                if (!alt_correct && base_correct) {
                    stats.altProvBaseWouldHaveCorrect++;
                }
            } else if (tage_provider == LONGEST) {

                if (LongestMatchPred == alttaken) {
                    stats.tageProvAltSame++;
                }
                if (LongestMatchPred == base_pred) {
                    stats.tageProvBaseSame++;
                }

                if (!long_correct && alt_correct) {
                    stats.tageProvAltWouldHaveCorrect++;
                }
                if (!long_correct && base_correct) {
                    stats.tageProvBaseWouldHaveCorrect++;
                }
            } else {

                if (base_pred == alttaken) {
                    stats.baseProvAltSame++;
                }
                if (base_pred == LongestMatchPred) {
                    stats.baseProvTageSame++;
                }

                if (!base_correct && long_correct) {
                    stats.baseProvTageWouldHaveCorrect++;
                }
                if (!base_correct && alt_correct) {
                    stats.baseProvAltWouldHaveCorrect++;
                }
            }
        } else {
            if (tage_provider == LONGEST) {
                if (LongestMatchPred == base_pred) {
                    stats.tageProvBaseSame++;
                }

                if (!tage_correct && base_correct) {
                    stats.tageProvBaseWouldHaveCorrect++;
                }

            } else {
                if (base_pred == LongestMatchPred) {
                    stats.baseProvTageSame++;
                }
                if (!base_correct && long_correct) {
                    stats.baseProvTageWouldHaveCorrect++;
                }
            }
        }
    }

    if (!tage_correct) stats.tageMispred++;
    if (!base_correct) stats.baseMispred++;
}


void TageBase::PrintStat(double instr) {

    printf("AllBr:%lu, CondBr:%i, TakenBr:%i, Ticks:%i\n",
            branchCount, stats.condBranches, stats.takenBranches,
            ticks
            );

    printf("Bi : %i(%.5f) | BiW : %i(%.5f) | BiC : %i(%.5f)\n",
            (stats.bimodalCorrect + stats.bimodalWrong),
           (double)(stats.bimodalCorrect + stats.bimodalWrong) / (double)stats.total,
           stats.bimodalWrong, (double)stats.bimodalWrong / (double)(stats.bimodalCorrect + stats.bimodalWrong),
           stats.bimodalCorrect, (double)stats.bimodalCorrect / (double)(stats.bimodalCorrect + stats.bimodalWrong)
            );
    printf("TG : %i(%.5f) | TGW : %i(%.5f) | TGC: %i(%.5f) \n",
            (stats.longestMatchCorrect + stats.longestMatchWrong),
           (double)(stats.longestMatchCorrect + stats.longestMatchWrong) / (double)stats.total,
           stats.longestMatchWrong, (double)stats.longestMatchWrong / (double)(stats.longestMatchCorrect + stats.longestMatchWrong),
           stats.longestMatchCorrect, (double)stats.longestMatchCorrect / (double)(stats.longestMatchCorrect + stats.longestMatchWrong));
    printf("ATG: %i(%.5f) | ATW : %i(%.5f) | ATC: %i(%.5f)\n",
            (stats.altMatchCorrect + stats.altMatchWrong),
           (double)(stats.altMatchCorrect + stats.altMatchWrong) / (double)stats.total,
           stats.altMatchWrong, (double)stats.altMatchWrong / (double)(stats.altMatchCorrect + stats.altMatchWrong),
           stats.altMatchCorrect, (double)stats.altMatchCorrect / (double)(stats.altMatchCorrect + stats.altMatchWrong));


    printf("Prov: T:[AS:%i AC:%i, BS:%i AC:%i], A:[TS:%i TC:%i, BS:%i BC:%i], B:[TS:%i TC:%i, AS:%i AC:%i]\n",
            stats.tageProvAltSame, stats.tageProvAltWouldHaveCorrect,
            stats.tageProvBaseSame, stats.tageProvBaseWouldHaveCorrect,
            stats.altProvTageSame, stats.altProvTageWouldHaveCorrect,
            stats.altProvBaseSame, stats.altProvBaseWouldHaveCorrect,
            stats.baseProvTageSame, stats.baseProvTageWouldHaveCorrect,
            stats.baseProvAltSame, stats.baseProvAltWouldHaveCorrect
            );
    printf("MPKI:: BASE:%.4f, TAGE:%.4f Red:%i (%.4f) Instr:%i\n",
            (double)stats.baseMispred / (double)instr * 1000,
            (double)stats.tageMispred / (double)instr * 1000,
            int(stats.baseMispred-stats.tageMispred),
            (double)(stats.baseMispred-stats.tageMispred) / (double)stats.baseMispred * 100,
            int(instr));


    for (int i = 0; i <= nhist; i++) {
        if (NOSKIP[i] || i == 0) {


        printf(
            "P%d: %i, %.5f | C: %i, %.5f | W: %i, %.5f | Alloc:%d | Allias: p:%i,n:%i\n",
            i, (stats.providerCorrect[i] + stats.providerWrong[i]),
            (double)(stats.providerCorrect[i] + stats.providerWrong[i]) / (double)stats.total,
            stats.providerCorrect[i],
            (double)stats.providerCorrect[i] / (double)stats.total,
            stats.providerWrong[i],
            (double)stats.providerWrong[i] / (double)stats.total,
            stats.allocations[i],
            stats.positiveAlliasing[i], stats.negativeAlliasing[i]
            );
        }
    }
    printf(
        "Allocations:%i, Init:%i, AperI:%.2f, TriesPerInit:%.3f | Useful: %i | TC: %i, TW: %i, uReset:%i\n",
        stats.totalAllocations,
        stats.totalAllocInit,
        (double)stats.totalAllocations / (double)stats.totalAllocInit,
        (double)stats.totalAllocTries / (double)stats.totalAllocInit,
        stats.totalUseful,
        (stats.longestMatchCorrect + stats.altMatchCorrect +
         stats.bimodalCorrect),
        (stats.longestMatchWrong + stats.altMatchWrong + stats.bimodalWrong),
        stats.uResets);

}




void TageBase::resetStats() {
    stats = {};

    for (int i = 0; i < (nbanklow * (1 << LogG)); i++) {
        Gentry &entry = gtable[1][i];
                entry.correct = 0;
                entry.incorrect = 0;
                entry.useful = 0;
    }
    for (int i = 0; i < (nbankhigh * (1 << LogG)); i++) {
        Gentry &entry = gtable[born][i];
                entry.correct = 0;
                entry.incorrect = 0;
                entry.useful = 0;
    }
};


}; // namespace LLBP