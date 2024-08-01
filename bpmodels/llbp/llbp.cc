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

#include "llbp.h"
#include "utils/intmath.hh"


namespace LLBP {


#define PRINTDEBUG 0
// #define COND (stats.total > 3000000)
#define COND (false)
#define COND2 (false)




#define HASH(...) HASHFN_(__VA_ARGS__)
#define HASHFN_(a1, a2, a3, a4) hashT(a1), hashN(a2), hashSk(a3), hashSh(a4)

#define TEST_ARGS 10, 1

#define CTX_RDIP

// #define INF_HIST
// #define MOVE_ON_ALLOC
// #define MOVE_USEFUL_ONLY
// #define MOVE_ON_USEFUL
// #define MOVE_ON_CORRECT
// #define MOVE_ON_EVICT
// #define MOVE_ON_USELESS

#ifdef LLBP_CONSTRAINED

#define L2_SORT_CTX
#define L2_SORT_PTRN

#define L2_SORT_CTX_NCONF
// #define L2_SORT_CTX_CORR

#endif


// #define L2_LRU
// #define PRINT_DISTANCE

#ifdef LLBP_CONSTRAINED
// #define LLBP_CONSTRAINED
// #define NO_ODD_TABLES
#define FILTER_TABLES
// #define NO_EVICT_LONGER
#else
// #define FILTER_TABLES

#endif

// #define CTX_BITS 14
// #define CTX_MASK uint64_t((1 << CTX_BITS) - 1)

// #define PRELOAD

// #define FILTER_TOPN_BRANCHES
// #define FILTER_H2P_FILTER
// #define USE_ORACLE

#define USE_NEW_REPLACE
// #define ALLOC_BOTH
// #define ALLOC_ONLY_L2
// #define USE_L2_FOR_TRAINING


#define L2_TABLE_START 0
#define L2_TABLE_END 40

// #define L2_TABLE_START 15
// #define L2_TABLE_END 30



// //================
// // As Victim cache
// #define USE_L2_AS_VICTIM
// #define MOVE_ON_EVICT

// #define DONT_MOVE_URESET
// #define DONT_MOVE_ALT_SAME
// #define DONT_MOVE_SIGN_CHANGE
// #define DONT_MOVE_BEFORE_USEFUL
#define DONT_MOVE_L2_USEFUL  //L1nRU
// #define REMOVE_USELESS_L2
// #define REMOVE_L2_IF_L1_SAME


#define OVERWRITE_SCL
// #define OVERWRITE_ALLWAYS



//================
// Alloc both
#define ALLOC_BOTH
#define ALLOC_SIMULTANEOUS
#define MOVE_ON_ALLOC

// //================
// // Use L2 for training
// #define USE_L2_FOR_TRAINING
// #define MOVE_ON_ALLOC

// //===============
// // Use only L2
// #define USE_ONLY_L2
// #define ALLOC_BOTH
// #define MOVE_ON_ALLOC


///================
// #define DONT_PREDICT_L2


/* Define the context hash function
 * 1. Type of history. Which branches should be hased
 *     0: All branches, 1: Only calls, 2: Calls and returns
 *     3: All unconditional branches, 4: All taken branches
 *
 * 2. Number of branches that should be hashed (W in the paper). Also, n
 *    or hashN in this code.
 * 3. Number of most recent branches skipped for CCID. Adds delay which
 *    is used to prefetch. Also, referred to as skip or HashSk in this code.
 * 4. Number of bits to shift the PC's. Is useful to avoid ping-pong context
 *    due to the XOR function in case a loop is executed
 *
 * =========================================================================
 * EXAMPLE                                                                 *
 *                                                                         *
 *                       pb-index (2.)  (3.)                               *
 *                      v             v     v                              *
 * history buffer : |l|k|j|i|h|g|f|e|d|c|b|a|                              *
 *                            ^prefetch (2.)^                              *
 *                                                                         *
 * a is the newest branch PC added to the buffer, l the oldest.            *
 * (2.) = n = W = 7; (3.) = 3                                              *
 * branches used to obtain PB index hash: j to d                           *
 * branches used to obtain hash to prefetch into PB: g to a                *
 * =========================================================================
*/

// LLBP default
#define HASHVALS 3, 8, 8, 2


///================



LLBP::LLBP(LLBPConfig cfg)
    : TageSCL(cfg.tsclConfig),
    CtrWidth(cfg.CtrWidth),
    ReplCtrWidth(cfg.ReplCtrWidth),
    CtxReplCtrWidth(cfg.CtxReplCtrWidth),

    llbpStorage(cfg.numContexts, cfg.numPatterns,
                 cfg.ctxAssoc, cfg.ptrnAssoc),
    numContexts(cfg.numContexts),
    numPatterns(cfg.numPatterns),
    // insertPossition(cfg.insertDistance),
    HASH(HASHVALS),
    // localCtxCacheSize(cfg.l2cacheSize),
    patternBuffer(cfg.l2cacheSize, cfg.l2cacheAssoc),
    // localCtxCache(cfg.l2cacheSize),
    TTWidth(cfg.TTWidth),
    CTWidth(cfg.CTWidth),
    simulateTiming(cfg.simulateTiming),
    warmup(false),
    accessDelay(cfg.accessDelay),

    // h2pFilter(512, 4),

    primMispLength(0,36,36),
    llbpMispLength(0,36,36),
    primProvLength(0,36,36),
    llbpProvLength(0,36,36),
    numHistPerContext(1,41,20),
    numUsefulHistPerContext(1,41,20)
{
    cfg.print();
    printf("LLBP:: CTX Hash:[%i,%i,%i,%i]\n",
            hashT, hashN, hashSk, hashSh);

    printf("CD: ");
    llbpStorage.printCfg();
    printf("PS: ");
    llbpStorage.insert(0,0)->patterns.printCfg();
    llbpStorage.erase(0);
    printf("PB: ");
    patternBuffer.printCfg();

    assert((!simulateTiming || (cfg.l2cacheSize >= hashSk))
            || "Pattern buffer hold at least all prefetches.");

    // Prefill the preload queue. Number of entries is equal to the
    // preload distance plus one.
    assert(hashSk >= 0);
    bbv[0].resize(window);
    bbv[1].resize(window);



    init_predictor();

    int mllbp[MAXNHIST];
    for (int i = 1; i <= nhist; i++) {
        mllbp[i] = (i%2) ? m[i] : m[i]+2;

        fghrT1[i] = new FoldedHistoryFast(ghr, mllbp[i], TTWidth);
        fghrT2[i] = new FoldedHistoryFast(ghr, mllbp[i], TTWidth - 1);
    }


#ifdef FILTER_TABLES
    // LLBP does not provide for all different history lenghts in
    // TAGE a prediction only for the following once. Note this
    // are not the actual length but the table index in TAGE
    auto l = {6,10,13,14,15,16,17,18,  19,20,22,24,26,28,32,36};

    int n = 0;
    for (auto i : l) {
        printf("%i=>%i:%i:%i ", i, n, n % cfg.ptrnAssoc, mllbp[i]);
        fltTables[i]=n++;

        // auto bucket = n / cfg.ptrnAssoc;
        // printf("%i=>%i:%i:%i ", i, n, bucket, mllbp[i]);
        // fltTables[i] = ( n << ceilLog2(cfg.ptrnAssoc) ) | bucket;
        // n++;
    }


    printf("\n");

#endif //FILTER_TABLES




}

LLBP::~LLBP() {

    for (int i = 1; i <= nhist; i++) {
        delete[] fghrT1[i];
        delete[] fghrT2[i];
    }
}




void LLBP::init_predictor() {
    printf("LLBP branch predictor\n");

}


uint64_t LLBP::calcHash(vector<uint64_t> &vec, int n) {
    uint64_t hash = 0;
    for (int i = 0; (i < n) && (i < vec.size()); i++) {
        hash ^= vec[vec.size() - i - 1];
    }
    return hash & ((1 << CTWidth) - 1);
}

/*
 * Given the {n} number of branches staring from vec[end-start]
 * to vec[end-start-n-1] we create the hash function by shifting
 * each PC by {shift} number if bits i.e.
 *
 *   000000000000|  PC  |    :vec[end-start]
 * ^ 0000000000|  PC  |00    :vec[end-start-1]
 * ^ 00000000|  PC  |0000    :vec[end-start-2]
 *           .                     .
 *           .                     .
 *           .                     .
 * ^ |  PC  |000000000000    :vec[end-start-n-1]
 * ----------------------
 *       final hash value
 * */
uint64_t LLBP::calcHash(std::list<uint64_t> &vec, int n, int start, int shift) {
    uint64_t hash = 0;
    if (vec.size() < (start + n)) {
        return 0;
    }
    uint64_t sh = 0;
    auto it = vec.begin();
    std::advance(it, start);
    for (; (it != vec.end()) && (n > 0); it++, n--) {
        uint64_t val = *it;

        // Shift the value
        hash ^= val << uint64_t(sh);

        sh += shift;
        if (sh >= CTWidth) {
            sh -= uint64_t(CTWidth);
        }
    }
    return hash & ((1 << CTWidth) - 1);
}

uint64_t LLBP::getCurContext(uint64_t pc) {
    return contexts[1].cur & ((1 << CTWidth) - 1);
}




/////////////////////////////////////////////////////////////////////////////////
// Main TAGE and chooser overwrites
/////////////////////////////////////////////////////////////////////////////////


bool LLBP::predict(uint64_t pc) {


    // 1. The base prediction
    TageSCL::basePredict(pc);
    tage_provider = BASE;
    bim_pred = base_pred;
    bimConf = baseConf;

    // 2. Make the LLBP prediction
    L2Predict(pc);

    // If there was a hit in level 2 mark it.
    // Also override the base prediction. The TAGE arbiter
    // will decide which prediction to use.
    // Note the TAGE chooser will use altConf in its index

#ifdef DONT_PREDICT_L2
    l2.isProvider = false;
#else
    baseConf = l2.hit ? l2.conf : bimConf;
    base_pred = l2.hit ? l2.pred : bim_pred;
    l2.isProvider = l2.hit;
#endif

    // 3. The TAGE prediction
    // Tage will call the chooseProvider function
    // Which arbitrates between the TAGE and L2 predictions.
    tagePredict(pc);

    DPRINTF("Prov: [TAGE:%i, L2:%i]\n", tage_pred, l2.isProvider);

    tage_scl_pred = tage_pred;
    scl_provider = tage_provider;

    // 3. SCL prediction
    SCLPredict(pc);

    // 4. Choose the correct prediction
    provider = scl_provider;

#ifdef OVERWRITE_SCL
    if (l2.isProvider) {
        provider = BASE;
        return l2.pred;
    }
#endif

    return tage_scl_pred;
}


unsigned LLBP::chooseProvider() {

    bool chooseL2 = l2.hit;

#ifdef DONT_PREDICT_L2
    chooseL2 = false;
#endif

#ifndef OVERWRITE_ALLWAYS
    // If LLBP hits we don't use LLBP if it has a longer history
    if (chooseL2 && (l2.histLength < HitBank)) {
        chooseL2 = false;
        l2.shorter = true;
    }
#endif

    if (simulateTiming && !warmup && !l2.prefetched) {
        chooseL2 = false;
    }


    if (chooseL2) {
        AltBank = 0;
        altConf = baseConf = l2.conf;
        alttaken = base_pred = l2.pred;
        l2.isProvider = true;
        return BASE;
    }

    // Clear the provider bit if instead the main TAGE is used.
    l2.isProvider = false;


    // If the longest is somehow certain use its prediction.
    if (tageConf != LowConf)
        return LONGEST;

    // Use on low confidence if the USE_ALT_ON_NA is negative
    if (use_alt_on_na[idxChooser()] < 0) {
        return LONGEST;
    }

    return (AltBank > 0) ? ALT : BASE;
}





// bool LLBP::isNotUseful(bool taken) {
// #if !defined(USE_NEW_REPLACE) || defined(DONT_PREDICT_L2)
//     return Base::isNotUseful(taken);
// #endif

//     // If there was no hit in level 2 use the default algorithm.
//     if (!l2.hit) {
//         return TageSCL::isNotUseful(taken);
//     }
//     return false;

//     // If both the longest and alternate predictions where correct
//     // we can possible free the longest entry to use it for other
//     // predictions.
//     if ((alttaken == taken) && (LongestMatchPred == taken)) {
//         // We only clear if the alternate prediction has a
//         // high confidence.
//         if (altConf == HighConf) {
//             if (AltBank > 0) {
//                 return true;
//             }
//         }
//     }

//     // We also clear if this is because the base predictor has
//     // overwritten the prediction and it was correct to do so.
//     // There is no need to handle this prediction in the main TAGE.
//     // if (!l2.override) {
//         // if ((l2.pred == taken) && (LongestMatchPred == taken)) {
//         //     l2stats.clearUbit++;
//         //     return true;
//         // }
//     // }

//     return false;
// }

// bool LLBP::isUseful(bool taken) {
// #if !defined(USE_NEW_REPLACE) || defined(DONT_PREDICT_L2)
//     return Base::isUseful(taken);
// #endif
//     if (!l2.hit) {
//         return TageSCL::isUseful(taken);
//     }
//     return false;
// }


bool LLBP::isNotUseful(bool taken) {

#if !defined(DONT_PREDICT_L2)
    if (l2.hit) {
        return false;
    }
#endif
    // If there was no hit in level 2 use the default algorithm.
    return TageSCL::isNotUseful(taken);
}

bool LLBP::isUseful(bool taken) {

#if !defined(DONT_PREDICT_L2)
    if (l2.hit) {
        return false;
    }
#endif
    // If there was no hit in level 2 use the default algorithm.
    return TageSCL::isUseful(taken);
}


void LLBP::updateL2Usefulness(bool taken) {

    if (!l2.hit) return;

    bool l1_correct = HitBank
                    ? (LongestMatchPred == taken)
                    : (bim_pred == taken);
    bool l2_correct = (l2.pred == taken);

    // In case both predictions are the same we can clear the u-bit
    // of the longer, more precise history.
    if (l1_correct && l2_correct) {

        // If both hit on the same table.
        if (l2.histLength == HitBank) {

#ifdef REMOVE_L2_IF_L1_SAME

            if (HitEntry->u < (1 << uwidth) - 1)
                HitEntry->u++;
            HitEntry->uselessReason = -1;
            l2stats.l1useful++;


            ctrupdate(L2HitEntry->replace, false, ReplCtrWidth);

            // if (L2HitEntry->replace == 0)
            {
                evictPattern(HitContext->key, L2HitEntry);
                HitContext->patterns.erase(L2HitEntry->key);
            }
#endif

        } else if (l2.histLength > HitBank) {

#ifdef REMOVE_USELESS_L2
               HitContext->patterns.erase(L2HitEntry->key);
#endif
        } else {
            // if (l2.conf == HighConf) {
                if (HitBank) {
#ifndef DONT_PREDICT_L2
#ifndef DONT_MOVE_L2_USEFUL
                    if (HitEntry->u > 0)
                        HitEntry->u--;
                    HitEntry->uselessReason = 4;
#endif
#endif
                }
            // }

#ifdef REMOVE_L2_IF_L1_SAME

            evictPattern(HitContext->key, L2HitEntry);
            HitContext->patterns.erase(L2HitEntry->key);
#endif

        }

    }
#ifdef DONT_PREDICT_L2
    bool provL2 = l2.hit;
#else
    bool provL2 = l2.isProvider;
#endif

    // If Level 1 was provider, it was correct and
    // level 2 was incorrect this prediction was useful.
    if (!provL2 && l1_correct && !l2_correct) {
        if (HitBank) {
#ifndef DONT_PREDICT_L2
            if (HitEntry->u < (1 << uwidth) - 1)
                HitEntry->u++;
#endif
        }
    }


    // Same for level 2. If it was the provider and was correct
    // but level 1 not, it was useful.
    if (l2.hit && !l2.shorter && l2_correct && !l1_correct) {

        ctrupdate(L2HitEntry->replace, true, ReplCtrWidth);

    }

}


void LLBP::updateTables(uint64_t pc, bool resolveDir, bool predDir) {


    DPRINTIF(COND,"%s nM:%i, TAGE:[d:%i, conf:%i, prov:%d HitBank:%d], BASE:[d:%i, conf:%i, %i]\n",
            (resolveDir != tage_pred) ? "Misp" : "Corr",
            stats.tageMispred,
            tage_pred, tageConf, tage_provider, HitBank, base_pred, baseConf, HitIdx);


    // 1. Table allocation --------------------

    bool ALLOC = false;
    if (l2.isProvider) {
        // If LLBP was provider we allocate if the prediction was wrong
        // and the history length is shorter than the maximum.
        ALLOC = (tage_pred != resolveDir) & (l2.histLength < nhist);

    } else {

        // If the prediction came from TAGE, it was wrong and the history
        // length is shorter than the maximum we allocate.
        ALLOC = (tage_pred != resolveDir) & (HitBank < nhist);

        // If LLBP was actually correct, it was longer than TAGE
        // but it was not choosen as provider we don't allocate.
        if (l2.hit && (l2.pred == resolveDir) && !l2.shorter) {
            ALLOC = false;
        }

        // This comes from the TAGE update function (alternative prediction)
        if (HitBank > 0) {
            if ((tageConf == LowConf) && (LongestMatchPred == resolveDir)) {
                ALLOC = false;
            }

            updateChooser(resolveDir);
        }

    }


    // Do the actual allocation
    // In case LLBP was the provider we overwrite the history length
    // of the TAGE prediction.
    auto tmp2 = HitBank;
    if (l2.isProvider) {
        HitBank = l2.histLength;
    }
    int nalloc = 0;
    if (ALLOC) {
        nalloc = 1+nnn;
        DPRINTF("Alloc:%i,%i, HL:%i, L2:[H:%i,S:%i,P:%i,D:%i]\n",
                stats.totalAllocInit, stats.total, HitBank, l2.hit, l2.shorter, l2.isProvider, l2.pred);
    }

    allocateTables(nalloc, pc, resolveDir);
    HitBank = tmp2;

    // 2. The LLBP + TAGE table updates
    // We only update either the TAGE or the LLBP tables.
    L2Update(pc, resolveDir, predDir);

    // 3. Finally the statistical corrector
    SCLUpdate(pc, resolveDir, predDir);
}


// Specialized update function of the primary predictor
//
bool LLBP::tageUpdateL1(uint64_t pc, bool resolveDir) {

    bool update_base = false;
    bool move_out = false;

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

#ifdef MOVE_ON_CORRECT
        move_out = LongestMatchPred == resolveDir;
#endif

        // Do the actual counter update
        ctrupdate(HitEntry->ctr, resolveDir, cwidth);
        // sign changes: no way it can have been useful
        if (HitEntry->ctr == (resolveDir ? 0 : -1)) {
            HitEntry->u = 0;
            move_out = false;
        }

        // If both the longest and alternate predictions where correct
        // we can possible free the longest entry to use it for other
        // predictions.
        // We clear this entry by clearing the useful bit.
        if (isNotUseful(resolveDir)) {
            if (HitEntry->u == 1) {
                HitEntry->u--;
                move_out = false;
            }
        }
        // If the longest hit was correct but the alternative prediction
        // was not promote this entry to be useful.
        if (isUseful(resolveDir)) {
            if (HitEntry->u < (1 << uwidth) - 1)
                HitEntry->u++;
            HitEntry->useful++;
            move_out = true;
// #if defined(ALLOC_BOTH)
//             auto ctx = llbpStorage.get(getCurContext(pc));
//             if (ctx) {
//                 auto ptrn = ctx->patterns.get(HitEntry->key);
//                 if (ptrn) {
//                     ptrn->useful++;
//                 }
//             }
// #endif
        }
        DPRINTIF(COND,"TageUpdate: idx:%d, ctr:%i,u:%i,T:%i\n",
                HitBank, HitEntry->ctr, HitEntry->u, HitEntry->tag);

        if (move_out) {
            // In case an entry became useful or provided a correct predicton
            // try to move it into the S1 predictor.
            // If successful the u-bit is cleared to make it available
            // for new allocations.

#ifndef MOVE_ON_ALLOC
#if defined(MOVE_ON_CORRECT) || defined(MOVE_ON_USEFUL)
            // S1MoveOut(HitBank, pc, resolveDir);

            if (L2Allocate(HitBank, pc, resolveDir, true)) {
                HitEntry->u = 0;
            }
#endif
#endif
        }

    } else {
        update_base = true;
    }

    // END TAGE UPDATE
    return update_base;
}











/////////////////////////////////////////////////////////////////////////////////
// LLBP PREDICTOR
/////////////////////////////////////////////////////////////////////////////////

void LLBP::L2Predict(uint64_t pc) {

    HitIdx = 0;

    // Calculate indices and tags
    // We need to do this explicity because we perform the prediction
    // before the TAGE prediction.
    calcIndicesAndTags(pc);

    L2HitEntry = nullptr;
    l2.histLength = 0;
    l2.isProvider = false;

    for (int i = 1; i <= nhist; i++) {

        if (!NOSKIP[i]) continue;



        uint64_t hash = 0;

        // We don't use all history lengths. Only 16
        auto _i = fltTables.contains(i) ? fltTables[i] : i;
        hash |= uint64_t(_i);
        // Table index (10 bits)


#define MAKE_KEY

#ifdef MAKE_KEY
        auto _key = 0;
        _key = pc;
        // _key ^= (pc >> (abs(logg[i] - i) + 1));

        // _key ^= tag1FHist[i]->value ^ (tag2FHist[i]->value << 1);
        _key ^= fghrT1[i]->value ^ (fghrT2[i]->value << 1);
        // _key ^= fghrT1[i]->value;
        // _key = _idx;


        // _key &= ((1 << (TTWidth + 4 * (i >= born))) - 1);
        _key &= ((1 << TTWidth) - 1);

        hash |= uint64_t(_key) << 10ULL;
#else

        // Index
        auto _idx = 0;
        _idx = pc ^ (pc >> (abs(logg[i] - i) + 1));
        _idx ^= indexFHist[i]->value;

        _idx &= ((1 << (logg[i])) - 1);

        auto _tag = 0;
        _tag = pc;
        _tag ^= tag1FHist[i]->value ^ (tag2FHist[i]->value << 1);

        // auto _key = _tag;

        // _tag &= ((1 << (TB[i])) - 1);
        _tag &= ((1 << TTWidth) - 1);


        // hash |= uint64_t(indexFHist[i]->value & ((1 << (logg[i])) - 1)) << 6ULL;
        // // PC rest of the bits
        // hash |= uint64_t(pc) << (16ULL);

        hash |= uint64_t(_tag) << 6ULL;
        // uint64_t ls = 20;
        // hash |= uint64_t(_idx) << 20ULL;

        // auto _key = _idx ^ (_tag << 2);


        // _key &= ((1 << TTWidth) - 1);

        // hash |= uint64_t(_key) << 10ULL;

#endif


        KEY[i] = hash;

        // KEY[i] = makeKey(pc, i);
    }



    // Get the current context (CCID)

    // Try full context first
    auto ctx_key = getCurContext(pc);
    PRINTIF(COND2,"%i L2Predict: %lx\n", branchCount, ctx_key);
    HitContext = llbpStorage.get(ctx_key);


    if (HitContext) {
        for (int i = nhist; i > 0; i--) {
            if (NOSKIP[i]) {
                // if (hitOnZeroTag && (GTAG[i] == 0)) {
                //     infHistories[KEY[i]].tag = 0;
                //     infHistories[KEY[i]].length = i;
                // }
                L2HitEntry = HitContext->patterns.get(KEY[i]);


                if (L2HitEntry) {
                    HitIdx = i;
                    hitVal[1] = L2HitEntry->ctr;
                    l2.histLength = i;
                    break;
                }
            }
        }
    }


    if (HitIdx > 0) {
        PRINTIF(COND,"S1Hit:%i,GI:%i,GT:%i,c:%i\n",
        L2HitEntry->length, GI[L2HitEntry->length], GTAG[L2HitEntry->length],
        L2HitEntry->ctr);

        // l2.l2provider = true;
    }

    l2.hit = (HitIdx > 0);
    l2.pVal = l2.hit ? hitVal[1] : BIM-2;
    l2.pred = l2.hit ? l2.pVal >= 0 : base_pred;
    l2.conf = compConf(l2.pVal, (HitIdx == 0) ? 2 : CtrWidth);
    l2.shorter = !l2.hit;

    if (simulateTiming) {
        L2PredHit = patternBuffer.get(ctx_key);
        if (L2PredHit) {
            L2PredHit->locked = false;
        }
        l2.prefetched = l2.hit && L2PredHit;
    } else {
        l2.prefetched = (ticks - lastMispredict) >= hashSk;
    }

}



// PREDICTOR UPDATE
void LLBP::L2Update(uint64_t pc, bool resolveDir, bool predDir) {


    // Update ----------------------------------------------------
    bool base_correct = (resolveDir == base_pred);
    bool l2_correct = (resolveDir == l2.pred);
    bool l1_correct = HitBank ? (LongestMatchPred == resolveDir) :
                      AltBank ? (alttaken == resolveDir) : (bim_pred == resolveDir);
    bool l2_useful = l2_correct && !l1_correct;

    bool updateBim = false;
    bool updateL1 = true;
    bool updateL2 = false;

#ifdef DONT_PREDICT_L2
    updateL2 = l2.hit;
    updateL1 = true;
#else
    updateL2 = l2.isProvider;
    // updateL2 = l2.hit;
    updateL1 = !l2.isProvider;
    // updateL1 = true;

#ifdef REMOVE_L2_IF_L1_SAME
    if (l2.isProvider && (HitBank == l2.histLength)) {
        updateL1 = true;
    }
#endif


#endif


    // In case level 2 provided the prediction
    // update the switching tables
    if (updateL2) {
        ctrupdate(L2HitEntry->ctr, resolveDir, CtrWidth);

#ifdef L2_SORT_CTX_CORR
        ctrupdate(HitContext->replace, l2_correct, CtxReplCtrWidth);
#endif



        DPRINTIF(COND,"S2Update: idx:%d, HL:%d ctr:%i,u:%i,T:%i\n",
                HitIdx, L2HitEntry->length, L2HitEntry->ctr, 0, L2HitEntry->tag);

#ifdef L2_SORT_CTX_NCONF
        // This function updates updates the context replacement counter
        // - If a pattern becomes confident (correct prediction)
        //   the replacement counter is increased
        // - If a pattern becomes low confident (incorrect prediction)
        //   the replacement counter is decreased
        if (L2HitEntry->ctr == (resolveDir ? 1 : -2)) {
        // else if (L2HitEntry->ctr == (resolveDir ? 2 : -3)) {
            // entry became medium confident
            ctrupdate(HitContext->replace, true, CtxReplCtrWidth);
        }
        else if (L2HitEntry->ctr == (resolveDir ? -1 : 0)) {
        // else if (L2HitEntry->ctr == (resolveDir ? -2 : 1)) {
            // entry became low confident
            ctrupdate(HitContext->replace, false, CtxReplCtrWidth);
        }
#endif


        if (!l2_correct && (l2.conf == LowConf)) {
            // if (HitBank > 0) {
            //     ctrupdate(HitEntry->ctr, resolveDir, CtrWidth);
            // } else {
                updateBim = true;
            // }
        }
    }

    if (updateL1) {
        // Update the first level TAGE
        updateBim = tageUpdateL1(pc, resolveDir);
    }

    if (updateBim) {

        // If the prediction was from the base predictor, update it.
        TageSCL::baseUpdate(pc, resolveDir, predDir);

        DPRINTIF(COND,"BIMUp: ctr:%i\n", BIM);
    }

    // Usefulness ------------------------------------------------
    if (l2.hit) {
        updateL2Usefulness(resolveDir);
    }

//     // Update the L2 predictor
    if (simulateTiming && L2PredHit) {
        if (updateL2) {
            L2PredHit->dirty = true;
        }
        if (l2.hit) {
            L2PredHit->used = true;
        }
        if (l2_useful) {
            L2PredHit->useful = true;
        }
    }

}

bool LLBP::L2Allocate(int histLen, uint64_t pc, bool taken, bool from_l1) {



#ifdef USE_L2_AS_VICTIM
    auto& victim = gtable[histLen][GI[histLen]];
    auto ctx_key = victim.ctx;
    auto k = victim.key;
    auto tag = victim.tag;
    auto idx = victim.idx;
    histLen = victim.hlen;
    pc = victim.pc;
    taken = victim.ctr >= 0;

    bool move = true;
    switch (victim.uselessReason) {
#ifdef DONT_MOVE_BEFORE_USEFUL
        case 0: move = false; break;
#endif
#ifdef DONT_MOVE_SIGN_CHANGE
        case 1: move = false; break;
#endif
#ifdef DONT_MOVE_ALT_SAME
        case 2: move = false; break;
#endif
#ifdef DONT_MOVE_URESET
        case 3: move = false; break;
#endif
#ifdef DONT_MOVE_L2_USEFUL
        case 4: move = false; break;
#endif
        default:
            // assert(false);
        break;
    }
    if (!move) return false;

    DPRINTIF(COND,"L1Evict:%i,GI:%i,T:%i,K:%#x,C:%#x\n", histLen, GI[histLen], victim.tag, victim.key, victim.ctx);
#endif



#ifdef FILTER_TABLES

    if (!fltTables.contains(histLen)) {
        return false;
    }
#endif //FILTER_TABLES


    // DPRINTIF(COND,"L2Alloc:%i,GI:%i,GT:%i\n", histLen, GI[histLen], GTAG[histLen]);


#ifndef USE_L2_AS_VICTIM


    // Create context key and pattern key
    auto ctx_key = getCurContext(pc);
    auto k = KEY[histLen];

    // Make a tag without hashing the PC
    // auto tag = gtag(pc, histLen);
    // auto idx = gindex(pc, histLen);
    auto tag = GTAG[histLen];
    auto idx = GI[histLen];
#endif


    auto ctx = allocateNewContext(pc, ctx_key);


    // Get the pattern
    auto ptrn = ctx->patterns.get(k);

    if (ptrn) {

        if (from_l1) {
            ptrn->ctr = HitEntry->ctr;
            ptrn->dir = HitEntry->ctr >= 0;
        }

        return true;
    }

    // Sorting before allocation to find the victim

#ifdef L2_SORT_PTRN
    ctx->sortPatters(k);
#endif

#ifdef LLBP_CONSTRAINED

    ptrn = ctx->patterns.getVictim(k);
    evictPattern(ctx_key, ptrn);
#endif

    ptrn = ctx->patterns.insert(k);

    ptrn->length = histLen;
    ptrn->tag = tag;
    ptrn->idx = idx;
    ptrn->useful = 0;
    ptrn->correct = 0;
    ptrn->key = k;
    ptrn->pc = pc;

    if (from_l1) {
        ptrn->ctr = HitEntry->ctr;
        ptrn->dir = HitEntry->ctr >= 0;
    } else {
        ptrn->ctr = taken ? 0 : -1;
        ptrn->dir = taken;
    }

    return true;
}



LLBP::Context* LLBP::allocateNewContext(uint64_t pc, uint64_t ctx_key) {

    Context* ctx = llbpStorage.get(ctx_key);

    if (!ctx) {

#ifdef L2_SORT_CTX
        llbpStorage.sortContexts(ctx_key);
#endif

        // Ensure the victim is not in the L2 predictor
        ctx = llbpStorage.getVictim(ctx_key);
        if (ctx) {
            for (auto& pt : ctx->patterns.getMap()) {
                evictPattern(ctx_key, &(pt.second->second), true);
            }
        }

        if (simulateTiming && ctx) {

            if (patternBuffer.exists(ctx->key)) {
                patternBuffer.erase(ctx->key);
            }

            // Also invalidate all entries in the preload queue
            // with this key.
            auto n = prefetchQueue.size();
            prefetchQueue.erase(std::remove_if(
                prefetchQueue.begin(), prefetchQueue.end(),
                [ctx_key](auto& e) { return e.key == ctx_key; }),
                prefetchQueue.end());
        }

        // Allocate a new context in the main buffer.
        ctx = llbpStorage.insert(ctx_key, pc);

        if (simulateTiming && ctx) {
            // Newly allocated.
            PBEntry entry(ctx_key);
            entry.valid = true;
            entry.dirty = true;
            entry.newlyAllocated = true;
            installInPB(entry, true);
        }

    }

    return ctx;
}


void LLBP::evictPattern(uint64_t ctx_key, HistEntry* ptrn, bool ctx_evict) {
    if (!ptrn) return;

#ifdef NO_EVICT_LONGER
    if (ptrn->length > histLen) {
        l2stats.noEvictLonger++;
        return false;
    }
#endif
    // auto _i = ctx->patterns.index(k);
    // printf("Evict: TN:%lx index:%lx key:%lx, \n", histLen, _i, k);
    // setIndex.insert(_i);


    auto& v_pattern = allPatterns[ctx_key][ptrn->key];
    v_pattern.length = ptrn->length;
    v_pattern.key = ptrn->key;
    v_pattern.idx = ptrn->idx;
    v_pattern.pc = ptrn->pc;

    v_pattern.correct += ptrn->correct;
    v_pattern.incorrect += ptrn->incorrect;
    v_pattern.useful += ptrn->useful;
    v_pattern.evicted += 1;
    if (ctx_evict)
        v_pattern.evicted_ctx += 1;

}



void LLBP::prefetch() {

    if (warmup) return;

    // Its impossible that we see two prefetches per cycle.
    // assert(prefetchQueue.empty() || (prefetchQueue.back().prefetchtime != ticks));

    // Perform the prefetching -----
    // Calculate the hash from the head of the history
    // auto ctx_key = calcHash(bbv[0], hashN, 0, hashSh) & CTX_MASK;
    auto ctx_key = contexts[2].cur;

    PRINTIF(COND2,"%i/%i Prefetch: %lx -> ", ticks, branchCount, ctx_key);



    // First check the preload queue if this entry is already enquened.
    auto it = std::find_if(
        prefetchQueue.begin(), prefetchQueue.end(),
                        [ctx_key](auto& e)
                        { return (e.key == ctx_key) && e.valid;});

    if (it != prefetchQueue.end()) {
        PRINTIF(COND2," Hit in prefetchQueue %lx", it->key);
        l2stats.l2PFHitInQueue++;
    }

    // Second check if its already cached.
    else if (patternBuffer.exists(ctx_key)) {
        // Copy the entry from the cache to the preload queue
        PRINTIF(COND2," Hit in pattern cache");
        l2stats.l2PFHitInCache++;
        patternBuffer.touch(ctx_key);
    }

    // Finally check if the context is available in the LLBP
    // and needs to be prefetched.
    else if (llbpStorage.exists(ctx_key)) {
        PRINTIF(COND2," Hit in CI -> prefetch");
        l2stats.l2PFHitInCI++;
        auto& pf_entry = prefetchQueue.emplace_back(ctx_key);
        pf_entry.valid = true;
        pf_entry.prefetchtime = ticks;
    } else {
        PRINTIF(COND2," Miss");
    }
    PRINTIF(COND2,"\n");

    // assert(prefetchQueue.size() <= hashSk + 1);
}

void LLBP::squashPrefetchQueue(bool btbMiss) {
    lastMispredict = ticks;
    if (btbMiss)
        l2stats.pfDroppedBTBMiss += prefetchQueue.size();
    else
        l2stats.pfDroppedMispredict += prefetchQueue.size();
    prefetchQueue.clear();
    if (btbMiss)
        prefetch();
}

void LLBP::tickPrefetchQueue() {

    // Tick should be called before the prediction is made.

    // Install prefeches if the prefetch delay has passed.
    if (!prefetchQueue.empty()) {
        auto& pf_entry = prefetchQueue.front();

        // If the prefetch delay has passed
        if (ticks - pf_entry.prefetchtime >= (accessDelay-1)) {

            PRINTIF(COND2," Install in Cache: %lx\n", pf_entry.key);
            pf_entry.locked = true;
            installInPB(pf_entry);
            prefetchQueue.pop_front();
        }
    }
}

void LLBP::installInPB(PBEntry &entry, bool bypass) {


    // First get the victim from the L2 predictor
    auto victim = patternBuffer.getVictim(entry.key);

    if (victim) {

        // If the entry is locked due to ungoing prefetch. Don't install in
        // PB but in LLBP right away.
        if (victim->locked && bypass) {
            l2stats.pfDroppedLocked++;
            l2stats.l2cacheDirtyEvict++;
            return;
        }

        if (victim->dirty) l2stats.l2cacheDirtyEvict++;
        else l2stats.l2cacheCleanEvict++;
        PRINTIF(COND2," Evict: %lx\n", victim->key);
    }

    // Copy the prefetched pattern set into the PB
    patternBuffer.insert(entry);
}




int LLBP::allocate(int Tidx, uint64_t pc, bool taken) {

    int alloc_res = -1;
#if defined(USE_L2_FOR_TRAINING) || defined(USE_ONLY_L2)
#ifdef FILTER_TOPN_BRANCHES
    if (!tracePCs.contains(pc)) {
        alloc_res = TageBase::allocate(Tidx, pc, taken);
    }
#endif
#else
    alloc_res = TageBase::allocate(Tidx, pc, taken);
#endif

    // Get the newly allocated entry and mark with the context and key
    if (alloc_res > 0) {
        auto& entry = gtable[Tidx][GI[Tidx]];
        entry.ctx = getCurContext(pc);
        entry.key = KEY[Tidx];
        DPRINTIF(COND,"L1Alloc:%i,GI:%i,T:%i,K:%#x,C:%#x\n", Tidx, GI[Tidx], entry.tag, entry.key, entry.ctx);

    }

#ifndef ALLOC_BOTH
    if (alloc_res > 0) {
        // Allocation successful
        return alloc_res;
    }
#endif

#ifdef ALLOC_SIMULTANEOUS
    if (alloc_res <= 0) {
        // Allocation not successful -> we also don't allocate in the L2
        return alloc_res;
    }
#endif


    // Allocation in the main TAGE not possible.
    // Allocate instead entries in the S1

#ifdef MOVE_ON_ALLOC
        if (L2Allocate(Tidx, pc, taken)) {

#ifndef DONT_PREDICT_L2
            stats.allocations[Tidx]++;
            stats.totalAllocations++;
            return 1;
#endif
        }
#endif


    return alloc_res;
}


// int LLBP::promote(HistEntry* l2entry) {
//     l2stats.l1promote++;

//     int idx = l2entry->length;
//     auto& entry = gtable[idx][l2entry->idx];
//     if (entry.tag == l2entry->tag) {
//         l2stats.l1promoteExits++;
//         return 0;
//     }

//     if (entry.u != 0) {
//         // DPRINTIF(COND,"NoAlloc:%i,GI:%i,GT:%i,pc%lx\n", idx, GI[idx], GTAG[idx],entry.pc);
//         l2stats.l1promoteFail++;
//         return -1;
//     }

// //     if (entry.pc != 0) {
// //         if (!overwriteNotUseful)
// //             return -1;

// //         DPRINTIF(COND,"AOverwrite:%i,GI:%i,GT:%i,pc%lx\n", idx, GI[idx], GTAG[idx],entry.pc);
// //         stats.overwriteAlloc[idx]++;
// //     }

// // #define OPTREMP
// // // the replacement is optimized with a single u bit: 0.2 %
// // #ifdef OPTREMP
// //     if (abs(2 * entry.ctr + 1) > 3) {
// //         if (entry.ctr > 0)
// //             entry.ctr--;
// //         else
// //             entry.ctr++;
// //         return 0;
// //     }
// // #endif

//     // evict(entry, idx);
//     // int c = checkReuse(entry, idx, pc,taken);

//     DPRINTIF(COND,"Alloc:%i,GI:%i,GT:%i\n",
//                   idx, l2entry->idx, l2entry->tag);

//     entry.tag = l2entry->tag;
//     entry.pc = l2entry->pc;
//     entry.ctx = HitContext->key;
//     entry.key = l2entry->key;
//     // entry.ctx = ;

//     entry.ctr = saturate(l2entry->ctr, cwidth);

//     // entry.u = 0;
//     stats.allocations[idx]++;
//     stats.totalAllocations++;
//     l2stats.l1promoteSuccess++;
//     entry.newalloc = true;
//     entry.allocTime = allocID[idx]++;
//     if (entry.correct < 0) entry.correct = 0;

//     return 1;
// }




void LLBP::updateStats(bool resolveDir, bool predDir, uint64_t pc) {

    TageSCL::updateStats(resolveDir, predDir, pc);


    // Check if storing the last history would had been useful.
    auto correct = predDir == resolveDir;

    auto llbp_correct = l2.isProvider && (resolveDir == l2.pred);

    bool prim_correct = (scl_provider == STC) ? (sc_pred == resolveDir) :
                        HitBank ? (LongestMatchPred == resolveDir) :
                        AltBank ? (alttaken == resolveDir) : (bim_pred == resolveDir);


    bool llbp_useful = llbp_correct && !prim_correct;


    if (l2.hit) {
        if (l2.isProvider) {
            if (llbp_correct) {
                if (prim_correct)
                    l2stats.l2OverrideSameCorr++;
                else
                    l2stats.l2OverrideGood++;
            } else {
                if (prim_correct)
                    l2stats.l2OverrideBad++;
                else
                    l2stats.l2OverrideSameWrong++;
            }
            if (L2HitEntry->pc != pc) {
                if (llbp_correct)
                    l2stats.ovrPosAlias++;
                else {
                    l2stats.ovrNegAlias++;
                }
            }
        } else {
            l2stats.l2NoOverride++;
        }
    }




    // Hits for contexts and patterns
    if (HitContext) {
        l2stats.l2CtxHit++;
        if (L2HitEntry) {
            l2stats.l2PtrnHit++;
        }
    }


    if (l2.isProvider) {
        l2stats.l2Prov++;
        llbpProvLength.insert(l2.histLength);
        if (llbp_correct) {
            l2stats.l2Correct++;
            L2HitEntry->correct++;
            HitContext->correct++;
            if (llbp_useful) {

                L2HitEntry->useful++;
                HitContext->useful++;
                if (L2HitEntry->useful == 1) {
                    HitContext->usefulPtrns++;
                }
            }


        } else {
            l2stats.l2Wrong++;
            llbpMispLength.insert(l2.histLength);
            L2HitEntry->incorrect++;
            HitContext->incorrect++;
        }
    } else {
        switch (provider) {
            case LONGEST:
            case ALT:
                {
                l2stats.tageProv++;
                auto l = (provider == LONGEST) ? HitBank : AltBank;
                primProvLength.insert(l);
                if (correct) l2stats.tageCorrect++;
                else {
                    l2stats.tageWrong++;
                    primMispLength.insert(l);
                }
                }
                break;

            case LOOP:
            case STC:
                l2stats.sclProv++;
                if (correct) l2stats.sclCorrect++;
                else l2stats.sclWrong++;
                break;

            case BASE:
                l2stats.baseProv++;
                if (correct) {
                    l2stats.baseCorrect++;
                } else {
                    l2stats.baseWrong++;
                }
                break;
            default:
                break;
        }
    }

    if (l2.hit && !l2.isProvider) {
        if (l2.shorter) {
            l2stats.l2notBecauseShorter++;
        }

        if (!l2.prefetched) {
            l2stats.l2notBecauseNotPrefetched++;
        }
    }
}



// Predictor update ----------------------------------------
void LLBP::UpdatePredictor(uint64_t PC, bool resolveDir,
                                  bool predDir, uint64_t branchTarget) {
    // Update our own base predictor

    branchCount++;
    stats.condBranches++;
    if (resolveDir)
        stats.takenBranches++;

    updateTables(PC, resolveDir, predDir);


    updateHistory(PC, resolveDir, OPTYPE_JMP_DIRECT_COND, branchTarget);
    updateStats(resolveDir, predDir, PC);

    bool do_prefetch = false;
    if (simulateTiming && (resolveDir != predDir)) {
        squashPrefetchQueue();
        do_prefetch = true;
    }

    // pcHistory->push(PC, resolveDir);
    // pcHistory->push((PC << 1 | resolveDir), resolveDir);

    do_prefetch |= updateRuntimeHash(PC, OPTYPE_JMP_DIRECT_COND, resolveDir);

    if (simulateTiming && do_prefetch) {
        prefetch();
    }
}


void LLBP::TrackOtherInst(uint64_t PC, OpType opType, bool taken,
                                 uint64_t branchTarget) {

    TageSCL::TrackOtherInst(PC, opType, taken, branchTarget);

    // pcHistory->push(PC, taken);
    // pcHistory->push((PC << 1 | taken), taken);

    auto do_prefetch = updateRuntimeHash(PC, opType, taken);
    if (simulateTiming && do_prefetch) {
        PRINTIF(COND2,"%i/%i Prefetch: %lx/t:%i from UpdateOther -> ", ticks, branchCount, PC, opType);
        prefetch();
    }
}



void LLBP::updateGHist(const bool bit) {
    TageSCL::updateGHist(bit);
    for (uint32_t i = 1; i <= nhist; ++i) {
        fghrT1[i]->update();
        fghrT2[i]->update();
    }
}


bool LLBP::updateRuntimeHash(uint64_t pc, OpType opType, bool taken) {

    // Hash of all branches
    auto isCall = (opType == OPTYPE_CALL_DIRECT_UNCOND)
               || (opType == OPTYPE_CALL_INDIRECT_UNCOND)
               || (opType == OPTYPE_CALL_DIRECT_COND);


    switch (hashT) {
        case 0: // All branches
        bbv[0].push_front(pc);
        bbv[1].push_front(branchCount);
        break;

        case 1: // Only calls
        if (isCall) {
            bbv[0].push_front(pc);
            bbv[1].push_front(branchCount);
        }
        break;

        case 2: // Only calls and returns
        if (isCall || (opType == OPTYPE_RET_UNCOND)) {
            bbv[0].push_front(pc);
            bbv[1].push_front(branchCount);
        }
        break;

        case 3: // Only unconditional branches
        if (opType != OPTYPE_JMP_DIRECT_COND) {
            bbv[0].push_front(pc);
            bbv[1].push_front(branchCount);
        }
        break;

        case 4: // All taken branches
        if (taken) {
            bbv[0].push_front(pc);
            bbv[1].push_front(branchCount);
        }
        break;
    }



    PRINTIF(COND,"UH:%llx, %i, %i\n", pc, opType, taken);
    // If the size has changed the hash has changed
    bool changed = false;
    if (bbv[0].size() > window) {
        changed = true;

        // Resize the history
        bbv[0].pop_back();
        bbv[1].pop_back();


        // The current context.
        updateContext(contexts[1], calcHash(bbv[0], hashN, hashSk, hashSh));
        // The short context.
        updateContext(contexts[3], calcHash(bbv[0], 2, hashSk, hashSh));

        // The prefetch context.
        updateContext(contexts[2], calcHash(bbv[0], hashN, 0, hashSh));

// #ifdef USE_L2_PRED_CACHE
//         // Advance the predictor queue
//         updateL2Predictor();
// #endif
    }
    return changed;
}

void LLBP::updateContext(CTX_t &ctx, uint64_t cur) {

    if (ctx.cur != cur) {
        ctx.prev = ctx.cur;
        ctx.cur = cur;
        ctx.lastSwitch = stats.total;
    }
}


void LLBP::tick() {
    TageSCL::tick();
    if (simulateTiming) {
        tickPrefetchQueue();
    }
}

void LLBP::btbMiss() {
    if (simulateTiming) {
        squashPrefetchQueue(true);
    }
}

void LLBP::setState(bool _warmup) {
    warmup = _warmup;
}






void LLBP::PrintStat(double instr) {

    TageSCL::PrintStat(instr);

    // Analyze the branch context
    numHistPerContext.reset();
    numUsefulHistPerContext.reset();

    int nPattern = 0, nUseful = 0, nUseful2 = 0;
    int nCtx = 0, nCtxUseful = 0;

    for (auto& ctx_pair : llbpStorage.getMap()) {

        auto& ctx = ctx_pair.second->second;

        int nuseful = 0;
        int pos = 0;
        for (auto& pt : ctx.patterns.getMap()) {
            if (pt.second->second.useful > 0) {
                nuseful++;
            }
        }

        int n = ctx.patterns.size();
        numHistPerContext.insert(n);
        numUsefulHistPerContext.insert(nuseful);
        nPattern += n;
        nUseful += nuseful;
        if (nuseful) {
            nCtxUseful++;
        }
    }




    printf("LLBP:: CtxHit:%i(%.4f), PtrnHit:%i(%.4f)\n",
            l2stats.l2CtxHit, l2stats.l2CtxHit / (double)stats.total,
            l2stats.l2PtrnHit, l2stats.l2PtrnHit / (double)stats.total
            );


    printf("PROVIDER::  BIM:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f] \n",
            l2stats.baseProv, (double)l2stats.baseProv / (double)stats.total,
            l2stats.baseCorrect, (double)l2stats.baseCorrect / (double)l2stats.baseProv,
            l2stats.baseWrong, (double)l2stats.baseWrong / (double)l2stats.baseProv,
            (double)l2stats.baseWrong / (double)instr * 1000
            );

    printf("PROVIDER:: TAGE:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f], \n",
            l2stats.tageProv, (double)l2stats.tageProv / (double)stats.total,
            l2stats.tageCorrect, (double)l2stats.tageCorrect / (double)l2stats.tageProv,
            l2stats.tageWrong, (double)l2stats.tageWrong / (double)l2stats.tageProv,
            (double)l2stats.tageWrong / (double)instr * 1000);

    printf("PROVIDER::  SCL:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f], \n",
            l2stats.sclProv, (double)l2stats.sclProv / (double)stats.total,
            l2stats.sclCorrect, (double)l2stats.sclCorrect / (double)l2stats.sclProv,
            l2stats.sclWrong, (double)l2stats.sclWrong / (double)l2stats.sclProv,
            (double)l2stats.sclWrong / (double)instr * 1000);

    printf("PROVIDER:: LLBP:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f], \n",
            l2stats.l2Prov, (double)l2stats.l2Prov / (double)stats.total,
            l2stats.l2Correct, (double)l2stats.l2Correct / (double)l2stats.l2Prov,
            l2stats.l2Wrong, (double)l2stats.l2Wrong / (double)l2stats.l2Prov,
            (double)l2stats.l2Wrong / (double)instr * 1000);


    printf("LLBP:: CtxHit:%i, PtrnHit:%i, Provider:%i(%.4f), NoProvider:[Shorter:%i(%.4f), NoPrefetch:%i(%.4f)]\n",
            l2stats.l2CtxHit, l2stats.l2PtrnHit, l2stats.l2Prov, (double)l2stats.l2Prov / (double)l2stats.l2PtrnHit,
            l2stats.l2notBecauseShorter, (double)l2stats.l2notBecauseShorter / (double)l2stats.l2PtrnHit,
            l2stats.l2notBecauseNotPrefetched, (double)l2stats.l2notBecauseNotPrefetched / (double)l2stats.l2PtrnHit
            );


    printf("LLBP:: PB Prefetch:[HitInPfq:%i, HitInPB:%i, HitInCI:%i], dropped[locked:%i, misp:%i, btbmiss:%i]\n",
            l2stats.l2PFHitInQueue, l2stats.l2PFHitInCache, l2stats.l2PFHitInCI, l2stats.pfDroppedLocked, l2stats.pfDroppedMispredict, l2stats.pfDroppedBTBMiss
            );

    auto tot_evicts = l2stats.l2cacheDirtyEvict + l2stats.l2cacheCleanEvict;
    printf("LLBP:: PB Evict:[Clean:%i(%.3f) Dirty:%i(%.3f)]\n",
            l2stats.l2cacheCleanEvict, (double)l2stats.l2cacheCleanEvict / (double)tot_evicts,
            l2stats.l2cacheDirtyEvict, (double)l2stats.l2cacheDirtyEvict / (double)tot_evicts
            );


    printf("LLBP:: LLBPHits:[NoOv:%i, SameCorr:%i, SameWrong:%i, GoodOv:%i, BadOv:%i] Alias:[P:%i(%.4f),N:%i(%.4f)]\n",
            l2stats.l2NoOverride, l2stats.l2OverrideSameCorr, l2stats.l2OverrideSameWrong, l2stats.l2OverrideGood, l2stats.l2OverrideBad,
            l2stats.ovrPosAlias, l2stats.ovrPosAlias / (double)l2stats.l2Prov,
            l2stats.ovrNegAlias, l2stats.ovrNegAlias / (double)l2stats.l2Prov
            );


    auto tot_pattern = (numPatterns * numContexts);


    nCtx = llbpStorage.getMap().size();
    printf(
        "LLBP:: Utilization: Patterns:[Total:%i,Alloc:%i(%.4f),Useful:%i(%.4f)], Ctx:[Total:%i,Alloc:%i(%.4f),Useful:%i(%.4f)]\n",

        tot_pattern, nPattern, nPattern / (double)tot_pattern,
        nUseful, nUseful / (double)tot_pattern,

        numContexts, nCtx, nCtx / (double)numContexts,
        nCtxUseful, nCtxUseful / (double)numContexts
        );

#define PRINTHIST

#ifdef PRINTHIST


    printf("Hist Histories per context\n");
    printf("%s\n", numHistPerContext.print(true,true).c_str());

    printf("Hist Useful histories per context\n");
    printf("%s\n", numUsefulHistPerContext.print(true,true).c_str());

    printf("Hist primary mispredict length (incorrect)\n");
    printf("%s\n", primMispLength.print(true,true).c_str());

    printf("Hist LLBP mispredict length (incorrect)\n");
    printf("%s\n", llbpMispLength.print(true,true).c_str());

    printf("Hist primary provider length\n");
    printf("%s\n", primProvLength.print(true,true).c_str());

    printf("Hist LLBP provider length\n");
    printf("%s\n", llbpProvLength.print(true,true).c_str());


#endif


}

void LLBP::resetStats() {
    TageSCL::resetStats();
    l2stats = {};

    primMispLength.reset();
    llbpMispLength.reset();
    primProvLength.reset();
    llbpProvLength.reset();
    numHistPerContext.reset();
    numUsefulHistPerContext.reset();
}


};  // namespace LLBP
