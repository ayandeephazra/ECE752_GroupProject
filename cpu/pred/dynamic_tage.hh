/*
 * Copyright (c) 2014 The University of Wisconsin
 *
 * Copyright (c) 2006 INRIA (Institut National de Recherche en
 * Informatique et en Automatique  / French National Research Institute
 * for Computer Science and Applied Mathematics)
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* @file
 * Implementation of a Dynamic TAGE branch predictor. TAGE is a global-history based
 * branch predictor. It features a PC-indexed bimodal predictor and N
 * partially tagged tables, indexed with a hash of the PC and the global
 * branch history. The different lengths of global branch history used to
 * index the partially tagged tables grow geometrically. A small path history
 * is also used in the hash.
 *
 * All TAGE tables are accessed in parallel, and the one using the longest
 * history that matches provides the prediction (some exceptions apply).
 * Entries are allocated in components using a longer history than the
 * one that predicted when the prediction is incorrect.
 */

#ifndef __CPU_PRED_TAGE_DYNAMIC_HH__
#define __CPU_PRED_TAGE_DYNAMIC_HH__

#include <vector>
#include "params/DynamicTAGE.hh"
#include "cpu/pred/tage_base.hh"


namespace gem5
{
  
  namespace branch_prediction
  {
    
    class DynamicTAGE : public TAGEBase {
      
    public:
      DynamicTAGE(const DynamicTAGEParams &p);
      void init() override;
      void buildTageTables() override;
      bool tagePredict(ThreadID tid, Addr branch_pc, bool cond_branch, BranchInfo* bi) override;
      
      
    protected:
      unsigned entries_per_tile; //Number of entries per tile is fixed during runtime.
      unsigned no_of_tiles ; //Number of tiles is fixed during runtime
      std::vector<int> current_configuration_vector;// = {5, 5, 5, 5, 5, 5, 5};
      std::vector<int> current_config;// = {0, 0, 0, 0, 0, 0, 0}; 
      TageEntry ***gtable;

      //	    const uint8_t MaxNumberOfTilesInOneTable;
      /*	    struct Dynamic_TAGE_tile{
		    TageEntry *tile_entries; //Pointer to a dynamically allocated array of TageEntries. The number of entries is determined by entries_per_tile value.
		    int32_t tile_ID; //Tile ID can change during runtime. The tileID is gotten from the configuration vector.
		    std::vector<bool> associated_history_length; //Since a tile can be associated to any history length during runtime, this vector mimics the interconnect switch that specifies to which history length this tile is associated to.
		    
		    }*/
      // Prediction Structures
      
      // Tage Entry struct from TAGEBase
      
    public:
      
      /**
       * Dynamically Changes the TAGE Table size according to the configuration vector
       */
      virtual void dynamicallyConfigTAGETable(std::vector<int> configuration_vector);
      //int8_t getCtr(int hitBank, int hitBankIndex) const;
      
      int calcTileIndex(int tableIndex){
	      return tableIndex / no_of_tiles;
      }
      
      int calcIndexWithinTile(int tableIndex){
	        return tableIndex % no_of_tiles;
      }
      void handleAllocAndUReset(bool alloc, bool taken, BranchInfo* bi, int nrand) override;
      void handleUReset() override;
      void handleTAGEUpdate(Addr branch_pc, bool taken, BranchInfo* bi) override;
     
};

} // namespace branch_prediction
} // namespace gem5 
#endif // __CPU_PRED_DYNAMIC_TAGE_HH__

