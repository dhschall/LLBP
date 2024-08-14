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

#include <assert.h>
#include "utils/intmath.hh"


namespace LLBP {


#define PRINTDEBUG 0
// #define COND (stats.total > 3000000)
#define COND (false)
#define COND2 (false)


#ifdef LLBP_CONSTRAINED
// Only sort if the constrained version is used
#define L2_SORT_CTX
#define L2_SORT_PTRN
#define FILTER_TABLES

#endif

#define OVERWRITE_SCL

/* The hash function is defined by 4 paramenters
 *
 * T: Type of history (T). Which branches should be hased
 *     0: All branches, 1: Only calls, 2: Calls and returns
 *     3: All unconditional branches, 4: All taken branches
 *
 * W: Number of branches that should be hashed (W in the paper).
 * D: Number of most recent branches skipped for CCID. Adds delay which
 *    is used to prefetch. (D in the paper.)
 * S: Number of bits to shift the PC's. Is useful to avoid ping-pong context
 *    due to the XOR function in case a loop is executed
 *
 * ********************************************************************* *
 * EXAMPLE                                                               *
 *                                                                       *
 *                       pb-index (2.)  (3.)                             *
 *                      v             v     v                            *
 * history buffer : |l|k|j|i|h|g|f|e|d|c|b|a|                            *
 *                            ^prefetch (2.)^                            *
 *                                                                       *
 * a is the newest branch PC added to the buffer, l the oldest.          *
 * (2.) = W = 7; (3.) = D = 3                                            *
 * branches used to obtain PB index hash: j to d                         *
 * branches used to obtain hash to prefetch into PB: g to a              *
 * ********************************************************************* *
 */

// LLBP default hash values
// [T, W, D, S]
#define HASHVALS 3, 8, 8, 2



LLBP::LLBP(LLBPConfig cfg)
    : TageSCL(cfg.tsclConfig),
    llbpStorage(cfg.numContexts, cfg.numPatterns,
                 cfg.ctxAssoc, cfg.ptrnAssoc),
    rcr(HASHVALS,cfg.CTWidth),
    patternBuffer(cfg.pbSize, cfg.pbAssoc),

    numContexts(cfg.numContexts),
    numPatterns(cfg.numPatterns),
    TTWidth(cfg.TTWidth),
    CtrWidth(cfg.CtrWidth),
    ReplCtrWidth(cfg.ReplCtrWidth),
    CtxReplCtrWidth(cfg.CtxReplCtrWidth),

    simulateTiming(cfg.simulateTiming),
    warmup(false),
    accessDelay(cfg.accessDelay),
    primMispLength(0,36,36),
    llbpMispLength(0,36,36),
    primProvLength(0,36,36),
    llbpProvLength(0,36,36),
    numHistPerContext(1,17,16),
    numUsefulHistPerContext(1,17,16)
{
    printf("LLBP branch predictor configs -------\n");
    cfg.print();

    printf("CD: ");
    llbpStorage.printCfg();
    printf("PS: ");
    llbpStorage.allocate(0,0)->patterns.printCfg();
    llbpStorage.erase(0);
    printf("PB: ");
    patternBuffer.printCfg();

    assert((!simulateTiming || (cfg.pbSize >= rcr.D))
            || "Pattern buffer hold at least all prefetches.");

    int mllbp[MAXNHIST];
    for (int i = 1; i <= nhist; i++) {
        mllbp[i] = (i%2) ? m[i] : m[i]+2;

        fghrT1[i] = new FoldedHistoryFast(ghr, mllbp[i], TTWidth);
        fghrT2[i] = new FoldedHistoryFast(ghr, mllbp[i], TTWidth - 1);
    }


#ifdef FILTER_TABLES
    // LLBP does not provide for all different history lenghts in
    // TAGE a prediction only for the following once which where
    // empirically determined. Note this
    // are not the actual length but the table indices in TAGE.
    auto l = {6,10,13,14,15,16,17,18,  19,20,22,24,26,28,32,36};

    int n = 0;
    for (auto i : l) {
        // To reduce the complexity of the multiplexer LLBP groups
        // always four consecutive history lenght in one bucket.
        // As the pattern sets are implemented a set associative
        // structure the lower bits determine the set=bucket.
        // The `fltTable`-map not only filters the history lengths
        // but also maps each length the correct pattern set index.
        // E.e. for the four way associativity the following function
        // ensures that history length 6,10,13,14 gets assign
        // 0,4,8,12 with the lowest two bits 0b00. Thus, the set will
        // be the same.
        auto bucket = n / cfg.ptrnAssoc;
        fltTables[i] = ((n%cfg.ptrnAssoc) << ceilLog2(cfg.ptrnAssoc) ) | bucket;
        printf("%i=>%i:%i:%i:%i ", i, n, bucket, fltTables[i], mllbp[i]);
        n++;
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
    llbpPredict(pc);

    // If there was a hit in level 2 mark it.
    // Also override the base prediction. The TAGE arbiter
    // will decide which prediction to use.
    // Note the TAGE chooser will use altConf in its index
    baseConf = llbp.hit ? llbp.conf : bimConf;
    base_pred = llbp.hit ? llbp.pred : bim_pred;
    llbp.isProvider = llbp.hit;

    // 3. The TAGE prediction
    // Tage will call the `chooseProvider` function
    // Which arbitrates between the TAGE and L2 predictions.
    tagePredict(pc);

    DPRINTF("Prov: [TAGE:%i, L2:%i]\n", tage_pred, llbp.isProvider);

    tage_scl_pred = tage_pred;
    scl_provider = tage_provider;

    // 3. SCL prediction
    SCLPredict(pc);

    // 4. Choose the correct prediction
    provider = scl_provider;

#ifdef OVERWRITE_SCL
    if (llbp.isProvider) {
        provider = BASE;
        return llbp.pred;
    }
#endif

    return tage_scl_pred;
}


unsigned LLBP::chooseProvider() {

    bool chooseL2 = llbp.hit;

    // If LLBP hits, we don't use LLBP if it has a longer history.
    if (chooseL2 && (llbp.histLength < HitBank)) {
        chooseL2 = false;
        llbp.shorter = true;
    }

    // Don't override if the prefetch is too late.
    if (simulateTiming && !warmup && !llbp.prefetched) {
        chooseL2 = false;
    }

    if (chooseL2) {
        AltBank = 0;
        altConf = baseConf = llbp.conf;
        alttaken = base_pred = llbp.pred;
        llbp.isProvider = true;
        return BASE;
    }

    // Clear the provider bit if instead the main TAGE is used.
    llbp.isProvider = false;


    // If the longest is somehow certain use its prediction.
    if (tageConf != LowConf)
        return LONGEST;

    // Use on low confidence if the USE_ALT_ON_NA is negative.
    if (use_alt_on_na[idxChooser()] < 0) {
        return LONGEST;
    }

    return (AltBank > 0) ? ALT : BASE;
}

bool LLBP::isNotUseful(bool taken) {

    if (llbp.hit) {
        return false;
    }
    // If there was no hit in level 2 use the default algorithm.
    return TageSCL::isNotUseful(taken);
}

bool LLBP::isUseful(bool taken) {

    if (llbp.hit) {
        return false;
    }
    // If there was no hit in LLBP use the default algorithm.
    return TageSCL::isUseful(taken);
}

bool LLBP::llbpCorrect(bool taken) {
    return llbp.hit && (taken == llbp.pred);
}

bool LLBP::primCorrect(bool taken) {
    return (scl_provider == STC) ? (sc_pred == taken) :
                         HitBank ? (LongestMatchPred == taken) :
                         AltBank ? (alttaken == taken) : (bim_pred == taken);
}

bool LLBP::tageCorrect(bool taken) {
    return  HitBank ? (LongestMatchPred == taken) :
            AltBank ? (alttaken == taken) : (bim_pred == taken);
}

bool LLBP::llbpUseful(bool taken) {
    return llbp.isProvider && (taken == llbp.pred) && !primCorrect(taken);
}

void LLBP::updateL2Usefulness(bool taken) {

    if (!llbp.hit) return;

    auto llbp_correct = llbpCorrect(taken);
    bool prim_correct = tageCorrect(taken);

    // If Level 1 was provider, it was correct and
    // level 2 was incorrect this prediction was useful.
    if (!llbp.isProvider && prim_correct && !llbp_correct) {
        if (HitBank) {
            if (HitEntry->u < (1 << uwidth) - 1)
                HitEntry->u++;
        }
    }
}


void LLBP::updateTables(uint64_t pc, bool resolveDir, bool predDir) {


    DPRINTIF(COND,"%s nM:%i, TAGE:[d:%i, conf:%i, prov:%d HitBank:%d], BASE:[d:%i, conf:%i]\n",
            (resolveDir != tage_pred) ? "Misp" : "Corr",
            stats.tageMispred,
            tage_pred, tageConf, tage_provider, HitBank, base_pred, baseConf);


    // 1. Table allocation --------------------
    bool ALLOC = false;
    if (llbp.isProvider) {
        // If LLBP was provider we allocate if the prediction was wrong
        // and the history length is shorter than the maximum.
        ALLOC = (tage_pred != resolveDir) & (llbp.histLength < nhist);

    } else {

        // If the prediction came from TAGE, it was wrong and the history
        // length is shorter than the maximum we allocate.
        ALLOC = (tage_pred != resolveDir) & (HitBank < nhist);

        // If LLBP was actually correct, it was longer than TAGE,
        // but it was not chosen as provider, then we don't allocate.
        if (llbp.hit && (llbp.pred == resolveDir) && !llbp.shorter) {
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
    // of the TAGE prediction. This forces the TAGE allocation
    // to start allocating with the history length of the LLBP.
    auto tmp2 = HitBank;
    if (llbp.isProvider) {
        HitBank = llbp.histLength;
    }
    int nalloc = 0;
    if (ALLOC) {
        nalloc = 1+nnn;
        DPRINTF("Alloc:%i,%i, HL:%i, L2:[H:%i,S:%i,P:%i,D:%i]\n",
                stats.totalAllocInit, stats.total, HitBank, llbp.hit, llbp.shorter, llbp.isProvider, llbp.pred);
    }

    allocateTables(nalloc, pc, resolveDir);
    HitBank = tmp2;

    // 2. The LLBP + TAGE table updates
    // We only update either the TAGE or the LLBP tables.
    llbpUpdate(pc, resolveDir, predDir);

    // 3. Finally, the statistical corrector.
    SCLUpdate(pc, resolveDir, predDir);
}



/////////////////////////////////////////////////////////////////////////////////
// LLBP PREDICTOR
/////////////////////////////////////////////////////////////////////////////////

void LLBP::llbpPredict(uint64_t pc) {

    // Calculate indices and tags
    // We need to do this explicity because we perform the prediction
    // before the TAGE prediction.
    calcIndicesAndTags(pc);

    llbpEntry = nullptr;
    llbp = {};

    for (int i = 1; i <= nhist; i++) {

        if (!NOSKIP[i]) continue;

        // We don't use all history lengths. Only 16
        // By using the lower bits for the table number we can
        // use it to manage the assocativity of the different history lengths.
        auto _i = fltTables.contains(i) ? fltTables[i] : i;
        // Table index (10 bits)


        auto _key =  pc;
        _key ^= fghrT1[i]->value ^ (fghrT2[i]->value << 1);
        // Mask the patterns bits
        _key &= ((1 << TTWidth) - 1);

        KEY[i] = uint64_t(_key) << 10ULL | uint64_t(_i);

    }



    // Get the current context (CCID)
    auto ctx_key = rcr.getCCID();
    PRINTIF(COND2,"%i L2Predict: %lx\n", branchCount, ctx_key);
    HitContext = llbpStorage.get(ctx_key);


    if (HitContext) {
        for (int i = nhist; i > 0; i--) {
            if (NOSKIP[i]) {
                llbpEntry = HitContext->patterns.get(KEY[i]);


                if (llbpEntry) {
                    llbp.hit = i;
                    llbp.pVal = llbpEntry->ctr;
                    llbp.pred = llbp.pVal >= 0;
                    llbp.conf = compConf(llbp.pVal, CtrWidth);
                    llbp.histLength = i;
                    break;
                }
            }
        }
    }


    if (llbp.hit) {
        PRINTIF(COND,"S1Hit:%i,GI:%i,GT:%i,c:%i\n",
        llbpEntry->length, GI[llbpEntry->length], GTAG[llbpEntry->length],
        llbpEntry->ctr);
    }

    // In case of timing simulation we check if the entry was already
    // prefetched.
    if (simulateTiming) {
        pbEntry = patternBuffer.get(ctx_key);
        if (pbEntry) {
            pbEntry->locked = false;
        }
        llbp.prefetched = llbp.hit && pbEntry;
    } else {
        llbp.prefetched = true;
    }

}



// PREDICTOR UPDATE
void LLBP::llbpUpdate(uint64_t pc, bool resolveDir, bool predDir) {

    // Update ----------------------------------------------------
    bool updateBim = false;
    bool updateLLBP = llbp.isProvider;
    bool updateTAGE = !llbp.isProvider;


    // Only the providing component is updated.
    // If the prediction came from LLBP its pattern gets updated.
    if (updateLLBP) {
        ctrupdate(llbpEntry->ctr, resolveDir, CtrWidth);

        // This function updates the context replacement counter
        // - If a pattern becomes confident (correct prediction)
        //   the replacement counter is increased
        // - If a pattern becomes low confident (incorrect prediction)
        //   the replacement counter is decreased
        if (llbpEntry->ctr == (resolveDir ? 1 : -2)) {
            // entry became medium confident
            ctrupdate(HitContext->replace, true, CtxReplCtrWidth);
        }
        else if (llbpEntry->ctr == (resolveDir ? -1 : 0)) {
            // entry became low confident
            ctrupdate(HitContext->replace, false, CtxReplCtrWidth);
        }

        // If the prediction wrong update also the BIM
        if (!llbpCorrect(resolveDir) && (llbp.conf == LowConf)) {
            updateBim = true;
        }
    }

    // If the prediction was from the TAGE predictor update it.
    if (updateTAGE) {
        updateBim = tageUpdate(pc, resolveDir);
    }

    // The base predictor is sometimes updated if the confidence of the
    // prediction is low.
    if (updateBim) {

        // If the prediction was from the base predictor, update it.
        TageSCL::baseUpdate(pc, resolveDir, predDir);

        DPRINTIF(COND,"BIMUp: ctr:%i\n", BIM);
    }

    // Usefulness ------------------------------------------------
    if (llbp.hit) {
        updateL2Usefulness(resolveDir);
    }


    // Update the pattern buffers statistics
    // and dirty bits.
    if (simulateTiming && pbEntry) {
        if (updateLLBP) {
            pbEntry->dirty = true;
        }
        if (llbp.hit) {
            pbEntry->used = true;
        }
        if (llbpUseful(resolveDir)) {
            pbEntry->useful = true;
        }
    }
}

bool LLBP::llbpAllocate(int histLen, uint64_t pc, bool taken) {

#ifdef FILTER_TABLES
    if (!fltTables.contains(histLen)) {
        return false;
    }
#endif //FILTER_TABLES


    // Create context key and pattern key
    auto ctx_key = rcr.getCCID();
    auto k = KEY[histLen];

    auto ctx = allocateNewContext(pc, ctx_key);

    // Check if the pattern already exists in LLBP.
    auto ptrn = ctx->patterns.get(k);
    if (ptrn) {
        return true;
    }

    // No pattern found. Allocate a new one.
    // Sorting before allocation to find the victim

#ifdef L2_SORT_PTRN
    ctx->sortPatters(k);
#endif

#ifdef LLBP_CONSTRAINED
    ptrn = ctx->patterns.getVictim(k);
#endif

    ptrn = ctx->patterns.insert(k);

    ptrn->length = histLen;
    ptrn->useful = 0;
    ptrn->correct = 0;
    ptrn->key = k;
    ptrn->pc = pc;

    ptrn->ctr = taken ? 0 : -1;
    ptrn->dir = taken;

    return true;
}



LLBP::Context* LLBP::allocateNewContext(uint64_t pc, uint64_t ctx_key) {

    Context* ctx = llbpStorage.get(ctx_key);

    // If the context does not exist we allocate a new one.
    if (!ctx) {

#ifdef L2_SORT_CTX
        llbpStorage.sortContexts(ctx_key);
#endif

        // Ensure the victim is not in the L2 predictor
        ctx = llbpStorage.getVictim(ctx_key);

        // If the victim context is still in pattern buffer
        // we need to remove it.
        if (simulateTiming && ctx) {

            if (patternBuffer.exists(ctx->key)) {
                patternBuffer.erase(ctx->key);
            }

            // Also invalidate all entries in the prefetch queue
            // with this key.
            auto n = prefetchQueue.size();
            prefetchQueue.erase(std::remove_if(
                prefetchQueue.begin(), prefetchQueue.end(),
                [ctx_key](auto& e) { return e.key == ctx_key; }),
                prefetchQueue.end());
        }

        // Allocate a new context in the LLBP storage.
        ctx = llbpStorage.allocate(ctx_key, pc);

        if (simulateTiming && ctx) {
            // Put the newly allocated entry into the PB.
            PBEntry entry(ctx_key);
            entry.valid = true;
            entry.dirty = true;
            entry.newlyAllocated = true;
            installInPB(entry, true);
        }

    }

    return ctx;
}



int LLBP::allocate(int Tidx, uint64_t pc, bool taken) {

    int alloc_res = TageBase::allocate(Tidx, pc, taken);

    // Get the newly allocated entry and mark with the context and key
    if (alloc_res > 0) {
        auto& entry = gtable[Tidx][GI[Tidx]];
        entry.ctx = rcr.getCCID();
        entry.key = KEY[Tidx];
        DPRINTIF(COND,"L1Alloc:%i,GI:%i,T:%i,K:%#x,C:%#x\n", Tidx, GI[Tidx], entry.tag, entry.key, entry.ctx);
    }

    if (alloc_res <= 0) {
        // Allocation not successful -> we also don't allocate in the LLBP
        return alloc_res;
    }

    // Try allocating in the LLBP.
    if (llbpAllocate(Tidx, pc, taken)) {
        stats.allocations[Tidx]++;
        stats.totalAllocations++;
        return 1;
    }
    return alloc_res;
}

// Prefetch functionality

void LLBP::prefetch() {

    if (warmup) return;

    // Perform the prefetching -----
    // Calculate the hash from the head of the history.
    auto ctx_key = rcr.getPCID();
    PRINTIF(COND2,"%i/%i Prefetch: %lx -> ", ticks, branchCount, ctx_key);



    // First check the preload queue if this entry is already enqueued.
    auto it = std::find_if(
        prefetchQueue.begin(), prefetchQueue.end(),
                        [ctx_key](auto& e)
                        { return (e.key == ctx_key) && e.valid;});

    if (it != prefetchQueue.end()) {
        PRINTIF(COND2," Hit in prefetchQueue %lx", it->key);
        llbpstats.l2PFHitInQueue++;
    }

    // Second check if its already cached.
    else if (patternBuffer.exists(ctx_key)) {
        // Copy the entry from the cache to the preload queue.
        PRINTIF(COND2," Hit in pattern cache");
        llbpstats.l2PFHitInCache++;
        patternBuffer.touch(ctx_key);
    }

    // Finally check if the context is available in the LLBP
    // and needs to be prefetched.
    else if (llbpStorage.exists(ctx_key)) {
        PRINTIF(COND2," Hit in CI -> prefetch");
        llbpstats.l2PFHitInCI++;
        auto& pf_entry = prefetchQueue.emplace_back(ctx_key);
        pf_entry.valid = true;
        pf_entry.prefetchtime = ticks;
    } else {
        PRINTIF(COND2," Miss");
    }
    PRINTIF(COND2,"\n");

}

void LLBP::squashPrefetchQueue(bool btbMiss) {
    if (btbMiss)
        llbpstats.pfDroppedBTBMiss += prefetchQueue.size();
    else
        llbpstats.pfDroppedMispredict += prefetchQueue.size();
    prefetchQueue.clear();

    // Once all prefetches are squashed we trigger prefetches
    // for an upcoming context.
    if (btbMiss)
        prefetch();
}

void LLBP::tickPrefetchQueue() {

    // Tick should be called before the prediction is made.

    // Install prefetches if the prefetch delay has passed.
    if (!prefetchQueue.empty()) {
        auto& pf_entry = prefetchQueue.front();

        // If the prefetch delay has passed
        if (ticks - pf_entry.prefetchtime >= (accessDelay)) {

            PRINTIF(COND2," Install in Cache: %lx\n", pf_entry.key);
            pf_entry.locked = true;
            installInPB(pf_entry);
            prefetchQueue.pop_front();
        }
    }
}

void LLBP::installInPB(PBEntry &entry, bool bypass) {


    // First get the victim from the LLBP predictor
    auto victim = patternBuffer.getVictim(entry.key);

    if (victim) {

        // If the entry is locked due to ongoing prefetch. Don't install in
        // PB but in LLBP right away.
        if (victim->locked && bypass) {
            llbpstats.pfDroppedLocked++;
            llbpstats.l2cacheDirtyEvict++;
            return;
        }

        if (victim->dirty) llbpstats.l2cacheDirtyEvict++;
        else llbpstats.l2cacheCleanEvict++;
        PRINTIF(COND2," Evict: %lx\n", victim->key);
    }

    // Copy the prefetched pattern set into the PB
    patternBuffer.insert(entry);
}


void LLBP::updateStats(bool resolveDir, bool predDir, uint64_t pc) {

    TageSCL::updateStats(resolveDir, predDir, pc);


    // Check if storing the last history would have been useful.
    auto correct = predDir == resolveDir;

    auto llbp_correct = llbp.isProvider && (resolveDir == llbp.pred);

    bool prim_correct = (scl_provider == STC) ? (sc_pred == resolveDir) :
                        HitBank ? (LongestMatchPred == resolveDir) :
                        AltBank ? (alttaken == resolveDir) : (bim_pred == resolveDir);


    bool llbp_useful = llbp_correct && !prim_correct;


    if (llbp.hit) {
        if (llbp.isProvider) {
            if (llbp_correct) {
                if (prim_correct)
                    llbpstats.l2OverrideSameCorr++;
                else
                    llbpstats.l2OverrideGood++;
            } else {
                if (prim_correct)
                    llbpstats.l2OverrideBad++;
                else
                    llbpstats.l2OverrideSameWrong++;
            }
            if (llbpEntry->pc != pc) {
                if (llbp_correct)
                    llbpstats.ovrPosAlias++;
                else {
                    llbpstats.ovrNegAlias++;
                }
            }
        } else {
            llbpstats.l2NoOverride++;
        }
    }




    // Hits for contexts and patterns
    if (HitContext) {
        llbpstats.l2CtxHit++;
        if (llbpEntry) {
            llbpstats.l2PtrnHit++;
        }
    }


    if (llbp.isProvider) {
        llbpstats.l2Prov++;
        llbpProvLength.insert(llbp.histLength);
        if (llbp_correct) {
            llbpstats.l2Correct++;
            llbpEntry->correct++;
            HitContext->correct++;
            if (llbp_useful) {

                llbpEntry->useful++;
                HitContext->useful++;
                if (llbpEntry->useful == 1) {
                    HitContext->usefulPtrns++;
                }
            }


        } else {
            llbpstats.l2Wrong++;
            llbpMispLength.insert(llbp.histLength);
            llbpEntry->incorrect++;
            HitContext->incorrect++;
        }
    } else {
        switch (provider) {
            case LONGEST:
            case ALT:
                {
                llbpstats.tageProv++;
                auto l = (provider == LONGEST) ? HitBank : AltBank;
                primProvLength.insert(l);
                if (correct) llbpstats.tageCorrect++;
                else {
                    llbpstats.tageWrong++;
                    primMispLength.insert(l);
                }
                }
                break;

            case LOOP:
            case STC:
                llbpstats.sclProv++;
                if (correct) llbpstats.sclCorrect++;
                else llbpstats.sclWrong++;
                break;

            case BASE:
                llbpstats.baseProv++;
                if (correct) {
                    llbpstats.baseCorrect++;
                } else {
                    llbpstats.baseWrong++;
                }
                break;
            default:
                break;
        }
    }

    if (llbp.hit && !llbp.isProvider) {
        if (llbp.shorter) {
            llbpstats.l2notBecauseShorter++;
        }

        if (!llbp.prefetched) {
            llbpstats.l2notBecauseNotPrefetched++;
        }
    }
}



// Predictor update ----------------------------------------
void LLBP::UpdatePredictor(uint64_t PC, bool resolveDir,
                                  bool predDir, uint64_t branchTarget) {

    // Update the TAGE and LLBP predictors via the
    // base class. This Will also update the histories and statistics.
    Tage::UpdatePredictor(PC, resolveDir, predDir, branchTarget);

    // Only thing left is the update of the prefetch queue and
    // the context hash.
    bool do_prefetch = false;
    if (simulateTiming && (resolveDir != predDir)) {
        squashPrefetchQueue();
        do_prefetch = true;
    }

    // In the default LLBP predictor there will be no update of the
    // runtime hash for conditional branches. However, this model
    // supports different types of histories.
    do_prefetch |= rcr.update(PC, OPTYPE_JMP_DIRECT_COND, resolveDir);
    if (simulateTiming && do_prefetch) {
        prefetch();
    }
}


void LLBP::TrackOtherInst(uint64_t PC, OpType opType, bool taken,
                                 uint64_t branchTarget) {

    TageSCL::TrackOtherInst(PC, opType, taken, branchTarget);

    auto do_prefetch = rcr.update(PC, opType, taken);
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


/************************************************************
 * RCR Functionality
 */

LLBP::RCR::RCR(int _T, int _W, int _D, int _shift, int _CTWidth)
    : CTWidth(_CTWidth), T(_T), W(_W), D(_D), S(_shift)
{
    bb[0].resize(maxwindow);
    bb[1].resize(maxwindow);
    ctxs = {0, 0};
    printf("\n\nRCR: context hash config: [T:%i, W:%i, D:%i, S:%i, CTWidth:%i]\n",
            T, W, D, S, CTWidth);
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
uint64_t LLBP::RCR::calcHash(std::list<uint64_t> &vec, int n, int start, int shift) {
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

uint64_t LLBP::RCR::getCCID() {
    return ctxs.ccid & ((1 << CTWidth) - 1);
}

uint64_t LLBP::RCR::getPCID() {
    return ctxs.pcid & ((1 << CTWidth) - 1);
}


bool LLBP::RCR::update(uint64_t pc, OpType opType, bool taken) {

    branchCount++;
    // Hash of all branches
    auto isCall = (opType == OPTYPE_CALL_DIRECT_UNCOND)
               || (opType == OPTYPE_CALL_INDIRECT_UNCOND)
               || (opType == OPTYPE_CALL_DIRECT_COND);


    switch (T) {
        case 0: // All branches
        bb[0].push_front(pc);
        bb[1].push_front(branchCount);
        break;

        case 1: // Only calls
        if (isCall) {
            bb[0].push_front(pc);
            bb[1].push_front(branchCount);
        }
        break;

        case 2: // Only calls and returns
        if (isCall || (opType == OPTYPE_RET_UNCOND)) {
            bb[0].push_front(pc);
            bb[1].push_front(branchCount);
        }
        break;

        case 3: // Only unconditional branches
        if (opType != OPTYPE_JMP_DIRECT_COND) {
            bb[0].push_front(pc);
            bb[1].push_front(branchCount);
        }
        break;

        case 4: // All taken branches
        if (taken) {
            bb[0].push_front(pc);
            bb[1].push_front(branchCount);
        }
        break;
    }



    PRINTIF(COND,"UH:%llx, %i, %i\n", pc, opType, taken);
    // If the size has changed the hash has changed
    bool changed = false;
    if (bb[0].size() > maxwindow) {
        changed = true;

        // Resize the history
        bb[0].pop_back();
        bb[1].pop_back();


        // The current context.
        ctxs.ccid = calcHash(bb[0], W, D, S);
        // The prefetch context.
        ctxs.pcid = calcHash(bb[0], W, 0, S);
    }
    return changed;
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

    int nPattern = 0, nUseful = 0;
    int nCtx = 0, nCtxUseful = 0;

    for (auto& ctx_pair : llbpStorage.getMap()) {

        auto& ctx = ctx_pair.second->second;

        int nuseful = 0;
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


    printf("LLBP branch predictor stats -------\n");

    printf("LLBP:: CtxHit:%i(%.4f), PtrnHit:%i(%.4f)\n",
           llbpstats.l2CtxHit, llbpstats.l2CtxHit / (double)stats.total,
           llbpstats.l2PtrnHit, llbpstats.l2PtrnHit / (double)stats.total
            );


    printf("PROVIDER::  BIM:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f] \n",
           llbpstats.baseProv, (double)llbpstats.baseProv / (double)stats.total,
           llbpstats.baseCorrect, (double)llbpstats.baseCorrect / (double)llbpstats.baseProv,
           llbpstats.baseWrong, (double)llbpstats.baseWrong / (double)llbpstats.baseProv,
           (double)llbpstats.baseWrong / (double)instr * 1000
            );

    printf("PROVIDER:: TAGE:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f], \n",
           llbpstats.tageProv, (double)llbpstats.tageProv / (double)stats.total,
           llbpstats.tageCorrect, (double)llbpstats.tageCorrect / (double)llbpstats.tageProv,
           llbpstats.tageWrong, (double)llbpstats.tageWrong / (double)llbpstats.tageProv,
           (double)llbpstats.tageWrong / (double)instr * 1000);

    printf("PROVIDER::  SCL:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f], \n",
           llbpstats.sclProv, (double)llbpstats.sclProv / (double)stats.total,
           llbpstats.sclCorrect, (double)llbpstats.sclCorrect / (double)llbpstats.sclProv,
           llbpstats.sclWrong, (double)llbpstats.sclWrong / (double)llbpstats.sclProv,
           (double)llbpstats.sclWrong / (double)instr * 1000);

    printf("PROVIDER:: LLBP:[P:%i(%.4f), C:%i(%.4f), W:%i(%.4f) MPKI:%.4f], \n",
           llbpstats.l2Prov, (double)llbpstats.l2Prov / (double)stats.total,
           llbpstats.l2Correct, (double)llbpstats.l2Correct / (double)llbpstats.l2Prov,
           llbpstats.l2Wrong, (double)llbpstats.l2Wrong / (double)llbpstats.l2Prov,
           (double)llbpstats.l2Wrong / (double)instr * 1000);


    printf("LLBP:: CtxHit:%i, PtrnHit:%i, Provider:%i(%.4f), NoProvider:[Shorter:%i(%.4f), NoPrefetch:%i(%.4f)]\n",
           llbpstats.l2CtxHit, llbpstats.l2PtrnHit, llbpstats.l2Prov, (double)llbpstats.l2Prov / (double)llbpstats.l2PtrnHit,
           llbpstats.l2notBecauseShorter, (double)llbpstats.l2notBecauseShorter / (double)llbpstats.l2PtrnHit,
           llbpstats.l2notBecauseNotPrefetched, (double)llbpstats.l2notBecauseNotPrefetched / (double)llbpstats.l2PtrnHit
            );


    printf("LLBP:: PB Prefetch:[HitInPfq:%i, HitInPB:%i, HitInCI:%i], dropped[locked:%i, misp:%i, btbmiss:%i]\n",
           llbpstats.l2PFHitInQueue, llbpstats.l2PFHitInCache, llbpstats.l2PFHitInCI, llbpstats.pfDroppedLocked, llbpstats.pfDroppedMispredict, llbpstats.pfDroppedBTBMiss
            );

    auto tot_evicts = llbpstats.l2cacheDirtyEvict + llbpstats.l2cacheCleanEvict;
    printf("LLBP:: PB Evict:[Clean:%i(%.3f) Dirty:%i(%.3f)]\n",
           llbpstats.l2cacheCleanEvict, (double)llbpstats.l2cacheCleanEvict / (double)tot_evicts,
           llbpstats.l2cacheDirtyEvict, (double)llbpstats.l2cacheDirtyEvict / (double)tot_evicts
            );


    printf("LLBP:: LLBPHits:[NoOv:%i, SameCorr:%i, SameWrong:%i, GoodOv:%i, BadOv:%i] Alias:[P:%i(%.4f),N:%i(%.4f)]\n",
           llbpstats.l2NoOverride, llbpstats.l2OverrideSameCorr, llbpstats.l2OverrideSameWrong, llbpstats.l2OverrideGood, llbpstats.l2OverrideBad,
           llbpstats.ovrPosAlias, llbpstats.ovrPosAlias / (double)llbpstats.l2Prov,
           llbpstats.ovrNegAlias, llbpstats.ovrNegAlias / (double)llbpstats.l2Prov
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
    llbpstats = {};

    primMispLength.reset();
    llbpMispLength.reset();
    primProvLength.reset();
    llbpProvLength.reset();
    numHistPerContext.reset();
    numUsefulHistPerContext.reset();
}


};  // namespace LLBP
