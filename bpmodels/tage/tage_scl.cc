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

#include "tage_scl.h"

#include "bpmodels/components/counters.h"
#include "utils/common.h"

using namespace std;

namespace LLBP {



//////////////////////////////////////////////////////////////////////////
// DEFINES

// --- Global --- //
// #define SC    // 8.2 % if TAGE alone
#define IMLI  // 0.2 %
#define LOCALH // 2.7 %
#define GLOBALH


#define LOOPPREDICTOR  // loop predictor enable

// --- SC --- //

// The statistical corrector components

#define PERCWIDTH 6  // Statistical corrector  counter width 5 -> 6 : 0.6 %


#ifdef TAGE8k
#define BWH

//The two BIAS tables in the SC component
//We play with confidence here
#define LOGBIAS 7

#define BIAS
#define BIASSK
#define BIASBANK

#else // TAGE64k

// The three BIAS tables in the SC component
// We play with the TAGE  confidence here, with the number of the hitting bank
#define BIAS
#define BIASSK
#define BIASBANK

#define LOGBIAS 8
#define LOCALS         // enable the 2nd local history
#define LOCALT         // enables the 3rd local history

#endif

#define INDBIAS                                                                 \
    (                                                                           \
        ((((PC ^ (PC >> 2)) << 1)                                               \
        ^ ((tageConf==LowConf) & (LongestMatchPred != alttaken)))  << 1)      \
        + pred_inter  \
    ) & ((1 << LOGBIAS) - 1)

#define INDBIASSK                                                              \
    (((((PC ^ (PC >> (LOGBIAS - 2))) << 1) ^ ((tageConf==HighConf))) << 1) + pred_inter) & \
        ((1 << LOGBIAS) - 1)

#define INDBIASBANK                                                      \
    (pred_inter + (((HitBank + 1) / 4) << 4) + ((tageConf==HighConf) << 1) +         \
     ((tageConf==LowConf) << 2) + ((AltBank != 0) << 3) + ((PC ^ (PC >> 2)) << 7)) & \
        ((1 << LOGBIAS) - 1)

// playing with putting more weights (x2)  on some of the SC components
// playing on using different update thresholds on SC
// update threshold for the statistical corrector
#define VARTHRES
#define WIDTHRES 12
#define WIDTHRESP 8
#ifdef VARTHRES
#define LOGSIZEUP 6  // not worth increasing
#else
#define LOGSIZEUP 0
#endif
#define LOGSIZEUPS (LOGSIZEUP / 2)
#define INDUPD (PC ^ (PC >> 2)) & ((1 << LOGSIZEUP) - 1)
#define INDUPDS ((PC ^ (PC >> 2)) & ((1 << (LOGSIZEUPS)) - 1))
#define EWIDTH 6



// --- Loop Predictor --- //

#ifdef LOOPPREDICTOR
// parameters of the loop predictor
#define WIDTHNBITERLOOP \
    10              // we predict only loops with less than 1K iterations
#define LOOPTAG 10  // tag width in the loop predictor

#endif

///////////////////////////////////////////////////////////////////////////

TageSCL::TageSCL(TSCLConfig cfg)
    : TageBase(cfg.tageConfig)
    , LogL(cfg.LogL)
    , useSC(cfg.useSC)
    , useLoop(cfg.useLoop)
    , disableConfCounter(false)
{
    cfg.print();
    init_predictor();
    predictorsize();
}

TageSCL::~TageSCL() {
    delete[] ltable;
}

void TageSCL::predictorsize() {
    int inter, STORAGESIZE = 0;

    if (useLoop) {

        inter = (1 << LogL) * (2 * WIDTHNBITERLOOP + LOOPTAG + 4 + 4 + 1);
        printf(" (LOOP %d) ", inter);
        STORAGESIZE += inter;
    }

    if (!useSC) {
        printf(" (SC %d) ", 0);
        return;
    }

    inter += WIDTHRES;
    inter = WIDTHRESP * ((1 << LOGSIZEUP));  // the update threshold counters
    inter +=
        3 * EWIDTH * (1 << LOGSIZEUPS);  // the extra weight of the partial sums
    inter += (PERCWIDTH) * 3 * (1 << (LOGBIAS));

    inter += (GNB - 2) * (1 << (LOGGNB)) * (PERCWIDTH) +
             (1 << (LOGGNB - 1)) * (2 * PERCWIDTH);
    inter += Gm[0];  // global histories for SC
    inter += (PNB - 2) * (1 << (LOGPNB)) * (PERCWIDTH) +
             (1 << (LOGPNB - 1)) * (2 * PERCWIDTH);
    // we use phist already counted for these tables

#ifdef BWH
    inter += BWNB * (1 << LOGBWNB) * PERCWIDTH;
    inter += EWIDTH * (1 << LOGSIZEUPS);	// the extra weight of the partial sums
    inter += BWm[0];
#endif

#ifdef LOCALH
    inter += (LNB - 2) * (1 << (LOGLNB)) * (PERCWIDTH) +
             (1 << (LOGLNB - 1)) * (2 * PERCWIDTH);
    inter += NLOCAL * Lm[0];
    inter += EWIDTH * (1 << LOGSIZEUPS);
#ifdef LOCALS
    inter += (SNB - 2) * (1 << (LOGSNB)) * (PERCWIDTH) +
             (1 << (LOGSNB - 1)) * (2 * PERCWIDTH);
    inter += NSECLOCAL * (Sm[0]);
    inter += EWIDTH * (1 << LOGSIZEUPS);

#endif
#ifdef LOCALT
    inter += (TNB - 2) * (1 << (LOGTNB)) * (PERCWIDTH) +
             (1 << (LOGTNB - 1)) * (2 * PERCWIDTH);
    inter += NTLOCAL * Tm[0];
    inter += EWIDTH * (1 << LOGSIZEUPS);
#endif

#endif

#ifdef IMLI

    inter += (1 << (LOGINB - 1)) * PERCWIDTH;
    inter += Im[0];

    inter += IMNB * (1 << (LOGIMNB - 1)) * PERCWIDTH;
    inter +=
        2 * EWIDTH * (1 << LOGSIZEUPS);  // the extra weight of the partial sums
    inter += 256 * IMm[0];
#endif
    inter += 2 * ConfWidth;  // the 2 counters in the choser
    STORAGESIZE += inter;

    printf(" (SC %d) ", inter);

}

void TageSCL::init_predictor() {

    ltable = new lentry[1 << (LogL)];


#ifdef LOOPPREDICTOR
    LVALID = false;
    WITHLOOP = -1;
    _lSeed = 0;
#endif

    updatethreshold = 35 << 3;

    for (int i = 0; i < (1 << LOGSIZEUP); i++) Pupdatethreshold[i] = 0;
    for (int i = 0; i < GNB; i++) GGEHL[i] = &GGEHLA[i][0];
    for (int i = 0; i < LNB; i++) LGEHL[i] = &LGEHLA[i][0];

    for (int i = 0; i < GNB; i++)
        for (int j = 0; j < ((1 << LOGGNB) - 1); j++) {
            if (!(j & 1)) {
                GGEHL[i][j] = -1;
            }
        }
    for (int i = 0; i < LNB; i++)
        for (int j = 0; j < ((1 << LOGLNB) - 1); j++) {
            if (!(j & 1)) {
                LGEHL[i][j] = -1;
            }
        }

    for (int i = 0; i < SNB; i++) SGEHL[i] = &SGEHLA[i][0];
    for (int i = 0; i < TNB; i++) TGEHL[i] = &TGEHLA[i][0];
    for (int i = 0; i < PNB; i++) PGEHL[i] = &PGEHLA[i][0];

    for (int i = 0; i < BWNB; i++) BWGEHL[i] = &BWGEHLA[i][0];
    for (int i = 0; i < BWNB; i++)
        for (int j = 0; j < ((1 << LOGBWNB) - 1); j++) {
            if (!(j & 1)) {
                BWGEHL[i][j] = -1;
            }
        }

#ifdef IMLI
#ifdef IMLIOH
    for (int i = 0; i < FNB; i++) FGEHL[i] = &FGEHLA[i][0];

    for (int i = 0; i < FNB; i++)
        for (int j = 0; j < ((1 << LOGFNB) - 1); j++) {
            if (!(j & 1)) {
                FGEHL[i][j] = -1;
            }
        }
#endif
    for (int i = 0; i < INB; i++) IGEHL[i] = &IGEHLA[i][0];
    for (int i = 0; i < INB; i++)
        for (int j = 0; j < ((1 << LOGINB) - 1); j++) {
            if (!(j & 1)) {
                IGEHL[i][j] = -1;
            }
        }
    for (int i = 0; i < IMNB; i++) IMGEHL[i] = &IMGEHLA[i][0];
    for (int i = 0; i < IMNB; i++)
        for (int j = 0; j < ((1 << LOGIMNB) - 1); j++) {
            if (!(j & 1)) {
                IMGEHL[i][j] = -1;
            }
        }

#endif
    for (int i = 0; i < SNB; i++)
        for (int j = 0; j < ((1 << LOGSNB) - 1); j++) {
            if (!(j & 1)) {
                SGEHL[i][j] = -1;
            }
        }
    for (int i = 0; i < TNB; i++)
        for (int j = 0; j < ((1 << LOGTNB) - 1); j++) {
            if (!(j & 1)) {
                TGEHL[i][j] = -1;
            }
        }
    for (int i = 0; i < PNB; i++)
        for (int j = 0; j < ((1 << LOGPNB) - 1); j++) {
            if (!(j & 1)) {
                PGEHL[i][j] = -1;
            }
        }

    for (int j = 0; j < (1 << LOGBIAS); j++) {
        switch (j & 3) {
            case 0:
                BiasSK[j] = -8;
                break;
            case 1:
                BiasSK[j] = 7;
                break;
            case 2:
                BiasSK[j] = -32;

                break;
            case 3:
                BiasSK[j] = 31;
                break;
        }
    }
    for (int j = 0; j < (1 << LOGBIAS); j++) {
        switch (j & 3) {
            case 0:
                Bias[j] = -32;

                break;
            case 1:
                Bias[j] = 31;
                break;
            case 2:
                Bias[j] = -1;
                break;
            case 3:
                Bias[j] = 0;
                break;
        }
    }
    for (int j = 0; j < (1 << LOGBIAS); j++) {
        switch (j & 3) {
            case 0:
                BiasBank[j] = -32;

                break;
            case 1:
                BiasBank[j] = 31;
                break;
            case 2:
                BiasBank[j] = -1;
                break;
            case 3:
                BiasBank[j] = 0;
                break;
        }
    }

    for (int i = 0; i < (1 << LOGSIZEUPS); i++) {
        WG[i] = 7;
        WL[i] = 7;
        WS[i] = 7;
        WT[i] = 7;
        WP[i] = 7;
        WI[i] = 7;
        WB[i] = 4;
	    WBW[i] = 7;
    }
    TICK = 0;
    for (int i = 0; i < NLOCAL; i++) {
        L_shist[i] = 0;
    }
    for (int i = 0; i < NSECLOCAL; i++) {
        S_slhist[i] = 0;
    }
}




void TageSCL::updateHistory(const uint64_t PC, const bool taken,
                              const OpType opType, const uint64_t branchTarget) {
    Tage::updateHistory(PC, taken, opType, branchTarget);

    if (opType != OPTYPE_JMP_DIRECT_COND) return;

#ifdef IMLI
    IMHIST[IMLIcount] = (IMHIST[IMLIcount] << 1) + taken;
    if (branchTarget < PC) {
        // This branch corresponds to a loop
        if (!taken) {
            // exit of the "loop"
            IMLIcount = 0;
        }
        if (taken) {
            if (IMLIcount < ((1 << Im[0]) - 1)) IMLIcount++;
        }
    }
#endif
    GHIST = (GHIST << 1) + (taken & (branchTarget < PC));
	BWHIST = (BWHIST << 1) + ((branchTarget < PC) & taken);
    L_shist[INDLOCAL] = (L_shist[INDLOCAL] << 1) + (taken);
    S_slhist[INDSLOCAL] = ((S_slhist[INDSLOCAL] << 1) + taken) ^ (PC & 15);
    T_slhist[INDTLOCAL] = (T_slhist[INDTLOCAL] << 1) + taken;
}

bool TageSCL::predict(uint64_t pc) {

    // 1. The base prediction
    basePredict(pc);
    tage_provider = BASE;

    // 2. The TAGE prediction
    tagePredict(pc);
    tage_scl_pred = tage_pred;
    scl_provider = tage_provider;

    // 3. SCL prediction
    SCLPredict(pc);

    // 4. Choose the correct prediction
    provider = scl_provider;
    return tage_scl_pred;
}

void TageSCL::SCLPredict(uint64_t pc) {

    // Loop prediction
    if (useLoop) {
        loop_pred = getloop(pc);  // loop prediction
        if ((WITHLOOP >= 0) && (LVALID)) {
            tage_scl_pred = loop_pred;
            scl_provider = LOOP;
        }
    }
    // Store the prediction without the SC
    pred_inter = tage_scl_pred;

    if (useSC) {
        int prov_inter = scl_provider;

        // Make the SC prediction
        SCpredict(pc);
        sc_pred = (LSUM >= 0);



        if (pred_inter != sc_pred) {
            // Chooser uses TAGE confidence and |LSUM|
            scl_provider = STC;

            // Minimal benefit in trying to avoid accuracy loss on low confidence SC
            // prediction and  high/medium confidence on TAGE
            //  but just uses 2 counters 0.3 % MPKI reduction
            if (!disableConfCounter){
                if ((tageConf==HighConf)) {
                    if ((abs(LSUM) < THRES / 4)) {
                        scl_provider = prov_inter;
                    }

                    else if ((abs(LSUM) < THRES / 2))
                        scl_provider = (SecondH < 0) ? STC : prov_inter;
                }

                if (tageConf==MedConf)
                    if ((abs(LSUM) < THRES / 4)) {
                        scl_provider = (FirstH < 0) ? STC : prov_inter;
                    }
            }
        }

        if (scl_provider == STC) {
            tage_scl_pred = sc_pred;
        }
    }
}

// PREDICTOR UPDATE
void TageSCL::updateTables(uint64_t PC, bool resolveDir, bool predDir)
{
    // TAGE update -----------------
    Tage::updateTables(PC, resolveDir, predDir);

    SCLUpdate(PC, resolveDir, predDir);
}

// PREDICTOR UPDATE
void TageSCL::SCLUpdate(uint64_t PC, bool resolveDir, bool predDir)
{
    // Loop update -----------------
    if (useLoop) {
        if (LVALID) {
            if (tage_scl_pred != loop_pred)
                ctrupdate(WITHLOOP, (loop_pred == resolveDir), 7);
        }
        loopupdate(PC, resolveDir, (tage_scl_pred != resolveDir));
    }

    // SC update -----------------
    if (useSC) {
        SCUpdate(PC, resolveDir, predDir);
    }
}



void TageSCL::SCpredict(uint64_t PC) {

    LSUM = 0;
    int8_t ctr = 0;

    LSUM += compPartial(PC);

    // integrate BIAS prediction -------
#ifdef BIAS
    ctr = Bias[INDBIAS];
    LSUM += (2 * ctr + 1);
#endif
#ifdef BIASSK
    ctr = BiasSK[INDBIASSK];
    LSUM += (2 * ctr + 1);
#endif
#ifdef BIASBANK
    ctr = BiasBank[INDBIASBANK];
    LSUM += (2 * ctr + 1);
#endif

    // Threshold for the statistical corrector
#ifdef VARTHRES
    LSUM = (1 + (WB[INDUPDS] >= 0)) * LSUM;
#endif

#ifdef GLOBALH
#ifdef TAGE8k
    LSUM += Gpredict (PC, GHIST, Gm, GGEHL, GNB, LOGGNB, WG);
    LSUM += Gpredict (PC, BWHIST, BWm, BWGEHL, BWNB, LOGBWNB, WBW);
#else // TAGE64k
    // integrate the GEHL predictions
    LSUM += Gpredict((PC << 1) + pred_inter, GHIST, Gm, GGEHL, GNB, LOGGNB, WG);
    LSUM += Gpredict(PC, phist, Pm, PGEHL, PNB, LOGPNB, WP);
#endif
#endif

    // Local history based components
#ifdef LOCALH
    LSUM += Gpredict(PC, L_shist[INDLOCAL], Lm, LGEHL, LNB, LOGLNB, WL);
#ifdef LOCALS
    LSUM += Gpredict(PC, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, LOGSNB, WS);
#endif
#ifdef LOCALT
    LSUM += Gpredict(PC, T_slhist[INDTLOCAL], Tm, TGEHL, TNB, LOGTNB, WT);
#endif
#endif

#ifdef IMLI
#ifndef TAGE8k
    LSUM += Gpredict(PC, IMHIST[(IMLIcount)], IMm, IMGEHL, IMNB, LOGIMNB, WIM);
#endif
    LSUM += Gpredict(PC, IMLIcount, Im, IGEHL, INB, LOGINB, WI);
#endif

    bool SCPRED = (LSUM >= 0);
    // Read the threshold for the statistical corrector
    // just  an heuristic if the respective contribution of component groups
    // can be multiplied by 2 or not
    THRES = (updatethreshold >> 3) + Pupdatethreshold[INDUPD]
#ifdef TAGE8k
#ifdef VARTHRES
      + 6 * ((WB[INDUPDS] >= 0)
#ifdef LOCALH
	     + (WL[INDUPDS] >= 0)
#endif
// #ifdef GSC
	     + (WG[INDUPDS] >= 0) + (WBW[INDUPDS] >= 0)
// #endif
#ifdef IMLI
	     + (WI[INDUPDS] >= 0)
#endif
      )
#endif

#else // TAGE64k
#ifdef VARTHRES
            +
            12 * ((WB[INDUPDS] >= 0) + (WP[INDUPDS] >= 0)
#ifdef LOCALH
                  + (WS[INDUPDS] >= 0) + (WT[INDUPDS] >= 0) + (WL[INDUPDS] >= 0)
#endif
                  + (WG[INDUPDS] >= 0)
#ifdef IMLI
                  + (WI[INDUPDS] >= 0)
#endif
                 )
#endif
#endif
        ;

}

void TageSCL::SCUpdate(uint64_t PC, bool resolveDir, bool predDir)
{

    bool SCPRED = (LSUM >= 0);

    // Chooser update -----------------
    if (!disableConfCounter) {
        if (pred_inter != SCPRED) {
            if ((abs(LSUM) < THRES))
                if (tageConf==HighConf) {
                    if ((abs(LSUM) < THRES / 2))
                        if ((abs(LSUM) >= THRES / 4))
                            ctrupdate(SecondH, (pred_inter == resolveDir),
                                    ConfWidth);
                }
            if (tageConf==MedConf)
                if ((abs(LSUM) < THRES / 4)) {
                    ctrupdate(FirstH, (pred_inter == resolveDir), ConfWidth);
                }
        }
    }


    if ((SCPRED != resolveDir) || ((abs(LSUM) < THRES))) {

        updatePartial(PC, resolveDir);

        {
            if (SCPRED != resolveDir) {
                Pupdatethreshold[INDUPD] += 1;
                updatethreshold += 1;
            }

            else {
                Pupdatethreshold[INDUPD] -= 1;
                updatethreshold -= 1;
            }

            if (Pupdatethreshold[INDUPD] >= (1 << (WIDTHRESP - 1)))
                Pupdatethreshold[INDUPD] = (1 << (WIDTHRESP - 1)) - 1;
            // Pupdatethreshold[INDUPD] could be negative
            if (Pupdatethreshold[INDUPD] < -(1 << (WIDTHRESP - 1)))
                Pupdatethreshold[INDUPD] = -(1 << (WIDTHRESP - 1));
            if (updatethreshold >= (1 << (WIDTHRES - 1)))
                updatethreshold = (1 << (WIDTHRES - 1)) - 1;
            // updatethreshold could be negative
            if (updatethreshold < -(1 << (WIDTHRES - 1)))
                updatethreshold = -(1 << (WIDTHRES - 1));
        }

#ifdef VARTHRES
        {
            int XSUM =
                LSUM - ((WB[INDUPDS] >= 0) *
                        ((2 * Bias[INDBIAS] + 1) + (2 * BiasSK[INDBIASSK] + 1) +
                         (2 * BiasBank[INDBIASBANK] + 1)));
            if ((XSUM + ((2 * Bias[INDBIAS] + 1) + (2 * BiasSK[INDBIASSK] + 1) +
                         (2 * BiasBank[INDBIASBANK] + 1)) >=
                 0) != (XSUM >= 0))
                ctrupdate(
                    WB[INDUPDS],
                    (((2 * Bias[INDBIAS] + 1) + (2 * BiasSK[INDBIASSK] + 1) +
                          (2 * BiasBank[INDBIASBANK] + 1) >=
                      0) == resolveDir),
                    EWIDTH);
        }
#endif
        // Bias update -----------------
        ctrupdate(Bias[INDBIAS], resolveDir, PERCWIDTH);
        ctrupdate(BiasSK[INDBIASSK], resolveDir, PERCWIDTH);
        ctrupdate(BiasBank[INDBIASBANK], resolveDir, PERCWIDTH);

        // Global history based components
#ifdef TAGE8k
        Gupdate(PC, resolveDir, GHIST, Gm, GGEHL, GNB, LOGGNB, WG);
        Gupdate(PC, resolveDir, BWHIST, BWm, BWGEHL, BWNB, LOGBWNB, WBW);
#else // TAGE64k
        Gupdate((PC << 1) + pred_inter, resolveDir, GHIST, Gm, GGEHL, GNB,
                LOGGNB, WG);
        Gupdate(PC, resolveDir, phist, Pm, PGEHL, PNB, LOGPNB, WP);
#endif

        // Local history based components
#ifdef LOCALH
        Gupdate(PC, resolveDir, L_shist[INDLOCAL], Lm, LGEHL, LNB, LOGLNB, WL);
#ifdef LOCALS
        Gupdate(PC, resolveDir, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, LOGSNB,
                WS);
#endif
#ifdef LOCALT

        Gupdate(PC, resolveDir, T_slhist[INDTLOCAL], Tm, TGEHL, TNB, LOGTNB,
                WT);
#endif
#endif

#ifdef IMLI
#ifndef TAGE8k
        Gupdate(PC, resolveDir, IMHIST[(IMLIcount)], IMm, IMGEHL, IMNB, LOGIMNB,
                WIM);
#endif
        Gupdate(PC, resolveDir, IMLIcount, Im, IGEHL, INB, LOGINB, WI);
#endif
    }
}





// #define GINDEX                                                                \
//     (((long long)PC) ^ bhist ^ (bhist >> (8 - i)) ^ (bhist >> (16 - 2 * i)) ^ \
//      (bhist >> (24 - 3 * i)) ^ (bhist >> (32 - 3 * i)) ^                      \
//      (bhist >> (40 - 4 * i))) &                                               \
//         ((1 << (logs - (i >= (NBR - 2)))) - 1)

int64_t TageSCL::gIndex(uint64_t PC, int64_t bhist, int logs, int nbr, int i)
{
    auto substr = (i >= (nbr - 2)) ? 1 : 0;
    return (((int64_t) PC) ^ bhist ^ (bhist >> (8 - i)) ^
                (bhist >> (16 - 2 * i)) ^ (bhist >> (24 - 3 * i)) ^
                (bhist >> (32 - 3 * i)) ^ (bhist >> (40 - 4 * i))) &
            ((1 << (logs - substr)) - 1);
}


int TageSCL::Gpredict(uint64_t PC, long long BHIST, int *length, int8_t **tab,
                        int NBR, int logs, int8_t *W) {
    int PERCSUM = 0;
    for (int i = 0; i < NBR; i++) {
        long long bhist = BHIST & ((long long)((1 << length[i]) - 1));
        long long index = gIndex(PC, bhist, logs, NBR, i);

        int8_t ctr = tab[i][index];

        PERCSUM += (2 * ctr + 1);
    }
#ifdef VARTHRES
    PERCSUM = (1 + (W[INDUPDS] >= 0)) * PERCSUM;
#endif
    return ((PERCSUM));
}

void TageSCL::Gupdate(uint64_t PC, bool taken, long long BHIST, int *length,
                        int8_t **tab, int NBR, int logs, int8_t *W) {
    int PERCSUM = 0;

    for (int i = 0; i < NBR; i++) {
        long long bhist = BHIST & ((long long)((1 << length[i]) - 1));
        long long index = gIndex(PC, bhist, logs, NBR, i);

        PERCSUM += (2 * tab[i][index] + 1);
        ctrupdate(tab[i][index], taken, PERCWIDTH);
    }
#ifdef VARTHRES
    {
        int XSUM = LSUM - ((W[INDUPDS] >= 0)) * PERCSUM;
        if ((XSUM + PERCSUM >= 0) != (XSUM >= 0))
            ctrupdate(W[INDUPDS], ((PERCSUM >= 0) == taken), EWIDTH);
    }
#endif
}

/// ---- LOOP PREDICTOR ---///

int TageSCL::lindex(uint64_t PC) {
    return (((PC ^ (PC >> 2)) & ((1 << (LogL - 2)) - 1)) << 2);
}

// loop prediction: only used if high confidence
// skewed associative 4-way
// At fetch time: speculative
#define CONFLOOP 15
bool TageSCL::getloop(uint64_t PC) {
    LHIT = -1;

    LI = lindex(PC);
    LIB = ((PC >> (LogL - 2)) & ((1 << (LogL - 2)) - 1));
    LTAG = (PC >> (LogL - 2)) & ((1 << 2 * LOOPTAG) - 1);
    LTAG ^= (LTAG >> LOOPTAG);
    LTAG = (LTAG & ((1 << LOOPTAG) - 1));

    for (int i = 0; i < 4; i++) {
        int index = (LI ^ ((LIB >> i) << 2)) + i;

        if (ltable[index].TAG == LTAG) {
            LHIT = i;
            LVALID = ((ltable[index].confid == CONFLOOP) ||
                      (ltable[index].confid * ltable[index].NbIter > 128));

            if (ltable[index].CurrentIter + 1 == ltable[index].NbIter)
                return (!(ltable[index].dir));
            return ((ltable[index].dir));
        }
    }

    LVALID = false;
    return (false);
}

void TageSCL::loopupdate(uint64_t PC, bool Taken, bool ALLOC) {
    if (LHIT >= 0) {
        int index = (LI ^ ((LIB >> LHIT) << 2)) + LHIT;
        // already a hit
        if (LVALID) {
            if (Taken != loop_pred) {
                // free the entry
                ltable[index].NbIter = 0;
                ltable[index].age = 0;
                ltable[index].confid = 0;
                ltable[index].CurrentIter = 0;
                return;

            } else if ((loop_pred != tage_pred) || ((MYLoopRANDOM() & 7) == 0))
                if (ltable[index].age < CONFLOOP) ltable[index].age++;
        }

        ltable[index].CurrentIter++;
        ltable[index].CurrentIter &= ((1 << WIDTHNBITERLOOP) - 1);
        // loop with more than 2** WIDTHNBITERLOOP iterations are not
        // treated correctly; but who cares :-)
        if (ltable[index].CurrentIter > ltable[index].NbIter) {
            ltable[index].confid = 0;
            ltable[index].NbIter = 0;
            // treat like the 1st encounter of the loop
        }
        if (Taken != ltable[index].dir) {
            if (ltable[index].CurrentIter == ltable[index].NbIter) {
                if (ltable[index].confid < CONFLOOP) ltable[index].confid++;
                if (ltable[index].NbIter < 3)
                // just do not predict when the loop count is 1 or 2
                {
                    // free the entry
                    ltable[index].dir = Taken;
                    ltable[index].NbIter = 0;
                    ltable[index].age = 0;
                    ltable[index].confid = 0;
                }
            } else {
                if (ltable[index].NbIter == 0) {
                    // first complete nest;
                    ltable[index].confid = 0;
                    ltable[index].NbIter = ltable[index].CurrentIter;
                } else {
                    // not the same number of iterations as last time: free
                    // the entry
                    ltable[index].NbIter = 0;
                    ltable[index].confid = 0;
                }
            }
            ltable[index].CurrentIter = 0;
        }

    } else if (ALLOC)

    {
        uint64_t X = MYLoopRANDOM() & 3;

        if ((MYLoopRANDOM() & 3) == 0)
            for (int i = 0; i < 4; i++) {
                int LHIT = (X + i) & 3;
                int index = (LI ^ ((LIB >> LHIT) << 2)) + LHIT;
                if (ltable[index].age == 0) {
                    ltable[index].dir = !Taken;
                    // most of mispredictions are on last iterations
                    ltable[index].TAG = LTAG;
                    ltable[index].NbIter = 0;
                    ltable[index].age = 7;
                    ltable[index].confid = 0;
                    ltable[index].CurrentIter = 0;
                    break;

                } else
                    ltable[index].age--;
                break;
            }
    }
}

void TageSCL::updateStats(bool taken, bool predtaken, uint64_t PC) {
    Tage::updateStats(taken, predtaken, PC);

    bool tage_correct = tage_pred == taken;
    bool tscl_correct = tage_scl_pred == taken;
    if (!tage_correct) sclstats.tageMisses++;
    if (!tscl_correct) sclstats.tsclMisses++;

    switch (provider) {
        case LONGEST:
        case ALT:
        case BASE:
            sclstats.provTage++;
            if (tage_scl_pred == taken) sclstats.tageCorrect++;
            else sclstats.tageIncorrect++;
            break;
        case LOOP:
            sclstats.provLoop++;
            if (tage_scl_pred == taken) sclstats.loopCorrect++;
            else sclstats.loopIncorrect++;
            break;
        case STC:
            sclstats.provSC++;
            if (tage_scl_pred == taken) sclstats.scCorrect++;
            else sclstats.scIncorrect++;
            break;
        default:
            break;
    }
}
void TageSCL::PrintStat(double instr) {
    Tage::PrintStat(instr);
    printf("SCL:: L:[P:%d(%.4f) C:%d(%.4f), W:%d(%.4f)] SC:[P:%d(%.4f) C:%d(%.4f), W:%d(%.4f)] TAGE:[P:%d(%.4f) C:%d(%.4f), W:%d(%.4f)] \n",
            sclstats.provLoop, sclstats.provLoop / (double)stats.total,
            sclstats.loopCorrect, sclstats.loopCorrect / (double)sclstats.provLoop,
            sclstats.loopIncorrect, sclstats.loopIncorrect / (double)sclstats.provLoop,
            sclstats.provSC,  sclstats.provSC / (double)stats.total,
            sclstats.scCorrect, sclstats.scCorrect / (double)sclstats.provSC,
            sclstats.scIncorrect, sclstats.scIncorrect / (double)sclstats.provSC,
            sclstats.provTage, sclstats.provTage / (double)stats.total,
            sclstats.tageCorrect, sclstats.tageCorrect / (double)sclstats.provTage,
            sclstats.tageIncorrect, sclstats.tageIncorrect / (double)sclstats.provTage);

    printf("MPKI:: TAGE:%.5f, SCL:%.5f Red:%.5f \n",
            (double)sclstats.tageMisses / (double)instr * 1000,
            (double)sclstats.tsclMisses / (double)instr * 1000,
            (double)(sclstats.tageMisses - sclstats.tsclMisses) / (double)sclstats.tageMisses * 100);
}

// Dump all tage tables as csv file

void TageSCL::DumpTables(std::string filename) {}


};  // namespace LLBP
