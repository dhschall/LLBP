/*
 * This file implements a basic Branch Target Buffer (BTB) structure, a Return Address Stack (RAS), and an indirect target branch prediction.
 * It was taken from the ChampSim framework provided with the paper https://doi.org/10.1109/IISWC59245.2023.00027.
 */

#include <iostream>

#include "basic_btb.h"
#include "ittage_64KB.h"
#include "ras.h"
#include "utils/common.h"
// #include <fmt/core.h>

#include <array>

class BTB
{

  // BasicBTB<8192, 8> btb[NUM_CPUS];
  BasicBTB<2048, 8> btb;
  my_predictor* ittage;
  RAS<64, 4096> ras;

  std::pair<uint64_t, uint8_t> interim_result;
  std::array<uint64_t, 12> btb_misses;


  public:

  void initialize()
  {
    std::cout << "BTB:" << std::endl;
    btb.initialize();
    std::cout << "Indirect:" << std::endl;
    ittage = new my_predictor();
    std::cout << "RAS:" << std::endl;
    ras.initialize();
    btb_misses.fill(0);
  }

  std::pair<uint64_t, uint8_t> prediction(uint64_t ip)
  {
    auto btb_pred = btb.predict(ip);
    interim_result = std::make_pair(0, false);
    if (btb_pred.first == 0 && btb_pred.second == BRANCH_INFO_ALWAYS_TAKEN) {
      // no prediction for this IP
    } else if (btb_pred.second == BRANCH_INFO_INDIRECT) {
      interim_result = std::make_pair(ittage->predict_brindirect(ip), true);
    } else if (btb_pred.second == BRANCH_INFO_RETURN) {
      interim_result = std::make_pair(ras.predict(), true);
    } else {
      interim_result = std::make_pair(btb_pred.first, btb_pred.second != BRANCH_INFO_CONDITIONAL);
    }
    return interim_result;
  }

  void update(uint64_t ip, uint64_t branch_target, uint8_t taken, OpType branch_type)
  {
    ittage->update_brindirect(ip, branch_type, taken, branch_target);
    ittage->fetch_history_update(ip, branch_type, taken, branch_target);
    ras.update(ip, branch_target, taken, branch_type);
    btb.update(ip, branch_target, taken, branch_type);

    if (taken && (branch_target != interim_result.first)) {
      btb_misses[branch_type]++;
    }
  }

  void final_stats(int inst)
  {

    std::array<std::pair<std::string, std::size_t>, 6> types{
        {std::pair{"BRANCH_DIRECT_JUMP", OPTYPE_JMP_DIRECT_UNCOND}, std::pair{"BRANCH_INDIRECT", OPTYPE_JMP_INDIRECT_UNCOND},
         std::pair{"BRANCH_CONDITIONAL", OPTYPE_JMP_DIRECT_COND}, std::pair{"BRANCH_DIRECT_CALL", OPTYPE_CALL_DIRECT_UNCOND},
         std::pair{"BRANCH_INDIRECT_CALL", OPTYPE_CALL_INDIRECT_UNCOND}, std::pair{"BRANCH_RETURN", OPTYPE_RET_UNCOND}}};

    // brpred->PrintStat(inst);
    // delete brpred;
    // printf("ZZZ Num_unique_taken_branches %lu\n", takenPCs.size());

    int total_misses = 0;
    for (auto [str, idx] : types) {
      printf("%s_BTB_MISS %i\n", str.c_str(), btb_misses[idx]);
      printf("%s_BTB_MPKI %.4f\n", str.c_str(), btb_misses[idx] * 1000.0 / (double)inst);
      total_misses += btb_misses[idx];
    }
    printf("TOTAL_BTB_MISS %i\n", total_misses);
    printf("TOTAL_BTB_MPKI %.3f\n", total_misses * 1000.0 / (double)inst);
  }
};