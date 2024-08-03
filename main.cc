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
 * Main file to simulate the LLBP branch predictor model.
 *
 * This framework aims to provide a fast and way to evaluate different branch
 * predictor configurations and explore the design space of LLBP.
 * It does not model the full pipeline but only the branch predictor. While
 * we approximate the timing impact of late prefetches by clocking the
 * predictor for every taken branch or if more than 8 instructions are
 * executed between branches this is only a rough estimation.
 * For the exact timing the predictor needs to be integrated with a full
 * CPU simulator like ChampSim or gem5.
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <boost/program_options.hpp>
#include <unordered_set>
using namespace std;

namespace po = boost::program_options;




#include "utils/common.h"
#include "utils/fileutils.h"
#include "bpmodels/base_predictor.h"
#include "btb/btb.cc"



#define COUNTER unsigned long long


#define d1K 1000
#define d10K 10000
#define d100K 100000
#define d1M 1000000
#define d10M 10000000
#define d30M 30000000
#define d60M 60000000
#define d100M 100000000
#define d300M 300000000
#define d600M 600000000
#define d1B 1000000000
#define d10B 10000000000
#define dMAX 0xFFFFFFFFFFFFFFFF

#define DEF_WARMUP 100*d1M
#define DEF_SIM 500*d1M

// #define WRITE_CSV
#define WRITE_CSV_START 5*d10M
#define WRITE_CSV_END  6*d10M
#define WRITE_CSV_COND_ONLY true

#define ANALYSIS



BasePredictor* brpred;


struct stats_t {
    uint64_t numCondMisp = 0;
    uint64_t numCond = 0;
    uint64_t numInstr = 0;
    uint64_t numTaken = 0;
    uint64_t numUncondBranches = 0;
    uint64_t numIndirect = 0;
};

stats_t total_stats;
stats_t roi_stats;
stats_t warmup_stats;

UINT64 nBranches = 0;
UINT64 instruction_count = 0;





std::string brmodel;
std::string outfile;
std::string trace_path;

bool tabledump = false;

uint64_t max_br_instruction = 0;
uint64_t warmup_instructions;
uint64_t sim_instructions;

// Filter never taken branches from BPU update
// This is the case for any commercial branch predictior.
// They update the BPU only for branches that are in the BTB which only
// contains taken branches.
std::unordered_set<uint64_t> takenPCs;
bool simulateBTB = false;


void CheckHeartBeat() {
    UINT64 dotInterval = 1000;
    UINT64 lineInterval = 1000 * dotInterval;

    if (nBranches % d1M == 0) {  // prints line every 1 million branches
        printf(
            "Nbr: %iM Ninst: %iM| MPKI: %10.4f : %llu,%llu | TotalMPKI: %10.4f : %llu,%llu | UC:%llu | nB:%llu\n",

            (int)(nBranches / d1M), (int)(instruction_count / d1M),
            1000.0 * (double)(total_stats.numCondMisp - warmup_stats.numCondMisp) / (double)(total_stats.numInstr - warmup_stats.numInstr),
            total_stats.numCondMisp - warmup_stats.numCondMisp, total_stats.numCond - warmup_stats.numCond,

            1000.0 * (double)total_stats.numCondMisp/ (double)total_stats.numInstr,
            total_stats.numCondMisp, total_stats.numCond,

            total_stats.numUncondBranches, nBranches
            );

        if (nBranches % (5*d1M) == 0) {
            printf("\n");
            brpred->PrintStat(instruction_count);
            printf("\n");
        }


        fflush(stdout);
    }

}  // void CheckHeartBeat


bool process_command_line(int argc, char** argv)
{

    try
    {
        po::options_description desc("Program Usage", 1024, 512);
        desc.add_options()
          ("help,h",         "produce help message")
          ("bpmodel,model,b",    po::value<std::string>(&brmodel), "define branch predictor model")
          ("input-file,i", po::value<std::string>(&trace_path)->required(), "input trace file")
          ("output,o",     po::value<std::string>(&outfile), "output file")
          ("simulate-btb",  po::bool_switch(&simulateBTB)->default_value(false), "Simulate BTB")
          ("tabledump,t",  po::bool_switch(&tabledump)->default_value(false), "dump TAGE tables")
          ("maxbrinst,m",  po::value<uint64_t>(&max_br_instruction)->default_value(dMAX), "max number of branches to simulate")
          ("inst-sim,n",   po::value<uint64_t>(&sim_instructions)->default_value(DEF_SIM), "max number of instructions to simulate")
          ("inst-warm,w",  po::value<uint64_t>(&warmup_instructions)->default_value(DEF_WARMUP), "number to instructions to warmup")
        ;

        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                    options(desc).positional(p).run(), vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }

        // There must be an easy way to handle the relationship between the
        // Yes, the magic is putting the po::notify after "help" option check
        po::notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    catch(...)
    {
        std::cerr << "Unknown error!" << "\n";
        return false;
    }

    return true;
}



int main(int argc, char* argv[]) {



    bool result = process_command_line(argc, argv);
    if (!result)
        return 1;


    // Extract the workload name from the path
    std::size_t pos2 = trace_path.rfind("/");
    std::size_t pos1 = trace_path.rfind('/', pos2 - 1);

    std::string workload = trace_path.substr(pos1 + 1, pos2 - pos1 - 1);

    printf("brmodel: %s, workload: %s %s trace: %s, output: %s\n",
            brmodel.c_str(), workload.c_str(), tabledump ? "tabledump " : "",
            trace_path.c_str(), outfile.c_str());

    ///////////////////////////////////////////////
    // Init variables
    ///////////////////////////////////////////////

    brpred = CreateBP(brmodel);

    ChampSimTrace trace;


    if (!trace.open(trace_path)) {
        printf("Error opening trace file: %s\n", trace_path.c_str());
        exit(-1);
    }

    BTB btb;
    if (simulateBTB) {
        btb.initialize();
    }


    ///////////////////////////////////////////////
    // read each trace record, simulate until done
    ///////////////////////////////////////////////

    OpType opType;
    UINT64 PC;
    bool branchTaken;
    UINT64 branchTarget;
    bool warmedUp = warmup_instructions ? false : true;

    branchTrace bt;

    // Skip the first branch/or header
    trace.getNextBranch(bt);

    printf("Starting simulation...\n");

    // while (trace.getNextBranchFromCSV(bt)) {
    while (trace.getNextBranch(bt)) {

        // printf("PC: %llx type: %x T %d N %d outcome: %d\n", bt.pc, (UINT32)bt.type, bt.taken, bt.skipped, bt.target);

        nBranches++;

        if (nBranches > max_br_instruction) {
            break;
        }

        if (instruction_count > (sim_instructions + warmup_instructions)) {
            break;
        }

        if (!warmedUp && (instruction_count >  warmup_instructions)) {
            warmedUp = true;
            warmup_stats = total_stats;
            printf("-----------------------------------------------\n"
                   "       Warmup complete after %llu instr \n"
                   "-----------------------------------------------\n",
                   warmup_instructions);
            brpred->resetStats();
        }



        try {
            opType = convertBrType(bt.type);
            PC = bt.pc;
            branchTaken = bt.taken;
            branchTarget = bt.target;
            instruction_count += bt.skipped + 1;
            total_stats.numInstr += bt.skipped + 1;

            CheckHeartBeat();

            branchTaken |= (opType != OPTYPE_JMP_DIRECT_COND);

            // In case we want to mimic the impact of late prefetches
            // we need to tick (clock) the predictor to simulate the prefetch
            // delay. As this framework only simulates the branch predictor
            // and not the full pipeline -- to speedup simulation and
            // exploration -- we don't know the exact number of cycles.
            // To get a rough estimation we tick the predictor for every
            // taken branch or if more than 8 instructions are skipped.
            // This is a rough estimation and based on the assumption that
            // the predictor can handle at most one taken branch per cycle,
            // an average instruction size of 4 bytes and a 32 byte
            // fetch buffer.
            // Note: For the paper we integrated the predictor with ChampSim
            // and we currently work on the integration with gem5.
            if (branchTaken || (bt.skipped > 8)) {
                brpred->tick();
            }

            /************************************************************************************************************/

            switch (opType) {
                case OPTYPE_OP:
                    break;

                // Conditional branches
                case OPTYPE_JMP_DIRECT_COND: {
                    bool predDir = false;

                    predDir = brpred->GetPrediction(PC);
                    brpred->UpdatePredictor(PC, branchTaken, predDir, branchTarget);

                    if (branchTaken) {
                        total_stats.numTaken++;
                    }

                    if (predDir != branchTaken) {
                        total_stats.numCondMisp++;
                    }

                    total_stats.numCond++;


                } break;

                case OPTYPE_JMP_INDIRECT_COND:
                case OPTYPE_JMP_INDIRECT_UNCOND:
                case OPTYPE_CALL_INDIRECT_UNCOND:
                    total_stats.numIndirect++;

                case OPTYPE_RET_UNCOND:
                case OPTYPE_CALL_DIRECT_UNCOND:
                case OPTYPE_JMP_DIRECT_UNCOND:
                    total_stats.numUncondBranches++;
                    // assert(branchTaken);
                    brpred->TrackOtherInst(PC, opType, branchTaken,
                                           branchTarget);
                    break;


                default:
                    fprintf(stderr, "OPTYPE_ERROR\n");
                    printf("OPTYPE_ERROR %i\n", opType);
                    exit(-1);  // this should never happen, if it does please
                               // email CBP org chair.
                    break;
            }

           if (simulateBTB) {
                // Simulate the BTB
                auto [predicted_branch_target, always_taken] = btb.prediction(PC);
                // Update the BTB
                btb.update(PC, branchTarget, branchTaken, opType);

                bool btbMiss = branchTaken && (branchTarget != predicted_branch_target);

                if (btbMiss) {
                    brpred->btbMiss();
                }
            }



            /************************************************************************************************************/
        } catch (const std::out_of_range& ex) {
            std::cout << ex.what() << '\n';
            break;
        }
    }

    ///////////////////////////////////////////
    // print_stats
    ///////////////////////////////////////////


    printf("  TRACE \t : %s\n", trace_path.c_str());
    printf("  TOTAL -------------------------------------\n");
    printf("  NUM_INSTRUCTIONS    : %10llu\n", instruction_count);
    printf("  NUM_BR              : %10llu (%.4f)\n", nBranches, (float)nBranches/(float)instruction_count);
    printf("  NUM_UNCOND_BR       : %10llu (%.4f)\n", total_stats.numUncondBranches, (float)total_stats.numUncondBranches/(float)nBranches);
    printf("  NUM_CONDITIONAL_BR  : %10llu (%.4f)\n", total_stats.numCond, (float)total_stats.numCond/(float)nBranches);
    printf("  NUM_TAKEN_BR        : %10llu (%.4f)\n", total_stats.numTaken, (float)total_stats.numTaken/(float)nBranches);
    printf("  NUM_MISPREDICTIONS  : %10llu\n", total_stats.numCondMisp);
    printf("  Cond Ratio          : %.6f\n",
           (double)(total_stats.numCondMisp) / (double)(nBranches));
    printf("  Cond MPKI           : %10.4f\n",
           1000.0 * (double)(total_stats.numCondMisp) / (double)(instruction_count));

    printf("  REGION OF INTEREST --------------------------\n");
    auto roi_instr = total_stats.numInstr - warmup_stats.numInstr;
    auto roi_misp = total_stats.numCondMisp - warmup_stats.numCondMisp;

    printf("  ROI INSTRUCTIONS   : %10llu\n", roi_instr);
    printf("  ROI MISPREDICTIONS : %10llu\n", roi_misp);
    printf("  ROI MPKI           : %10.4f\n", 1000.0 * (double)(roi_misp) / (double)(roi_instr));
    printf("\n");

    brpred->PrintStat(roi_instr);

    btb.final_stats(roi_instr);
    if (tabledump)
        brpred->DumpTables(outfile + "-table_dump.csv");

    printf("\n");

    delete brpred;
}
