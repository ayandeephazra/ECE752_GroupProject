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
 * Implementation of a TAGE branch predictor
 */

#include "cpu/pred/dynamic_tage.hh"

#include "base/intmath.hh"
#include "base/logging.hh"
#include "debug/Fetch.hh"
#include "debug/Tage.hh"
#include <numeric>


namespace gem5
{

  namespace branch_prediction
  {

    DynamicTAGE::DynamicTAGE(const DynamicTAGEParams &p)
      : TAGEBase(p){
		entries_per_tile = p.entries_per_tile;
		no_of_tiles = p.no_of_tiles;
		//int n = sizeof(p.current_config) / sizeof(p.current_config[0]);
		current_config = p.current_config;
		current_configuration_vector = p.current_configuration_vector;

		//fill(current_configuration_vector.begin(), current_configuration_vector.end(), 0);
	  	//fill(current_config.begin(), current_config.end(), 5);
    }

    // TAGEBase::BranchInfo*
    // TAGEBase::makeBranchInfo() {
    //     return new BranchInfo(*this);
    // }

    void
    DynamicTAGE::init()
    {
      //  TAGEBase::init();
      if (initialized) {
	return;
      }
    
      // Current method for periodically resetting the u counter bits only
      // works for 1 or 2 bits
      // Also make sure that it is not 0
      assert(tagTableUBits <= 2 && (tagTableUBits > 0));

      // we use int type for the path history, so it cannot be more than
      // its size
      assert(pathHistBits <= (sizeof(int)*8));

      // initialize the counter to half of the period
      assert(logUResetPeriod != 0);
      tCounter = initialTCounterValue;

      assert(histBufferSize > maxHist * 2);

      useAltPredForNewlyAllocated.resize(numUseAltOnNa, 0);

      for (auto& history : threadHistory) {
        history.pathHist = 0;
        history.globalHistory = new uint8_t[histBufferSize];
        history.gHist = history.globalHistory;
        memset(history.gHist, 0, histBufferSize);
        history.ptGhist = 0;
      }

      histLengths = new int [nHistoryTables+1];

      calculateParameters();

      //assert(tagTableTagWidths.size() == (nHistoryTables+1));
      //assert(logTagTableSizes.size() == (nHistoryTables+1));

      // First entry is for the Bimodal table and it is untagged in this
      // implementation
      assert(tagTableTagWidths[0] == 0);

      for (auto& history : threadHistory) {
        history.computeIndices = new FoldedHistory[nHistoryTables+1];
        history.computeTags[0] = new FoldedHistory[nHistoryTables+1];
        history.computeTags[1] = new FoldedHistory[nHistoryTables+1];

        initFoldedHistories(history);
      }

      const uint64_t bimodalTableSize = 1ULL << logTagTableSizes[0];
      btablePrediction.resize(bimodalTableSize, false);
      btableHysteresis.resize(bimodalTableSize >> logRatioBiModalHystEntries, true);

      gtable = new TageEntry** [nHistoryTables + 1];
      buildTageTables();
	  //current_config = 
	  //fill(current_configuration_vector.begin(), current_configuration_vector.end(), 0);
	  //fill(current_config.begin(), current_config.end(), 5);
	  dynamicallyConfigTAGETable(current_configuration_vector);


      tableIndices = new int [nHistoryTables+1];
      tableTags = new int [nHistoryTables+1];
      initialized = true;
    }

    void
    DynamicTAGE::buildTageTables() {
      /*
	here, we are now trying to allocate a bunch of array of pointers, to an array of pointers. 
	Thus, a 2D array of pointers is gotten. The pointers in each column correspond to a dynamic TAGE tile. 
	While, each column represents a dynamic TAGE table. */
      for (int i = 1; i <= nHistoryTables; i++) {
		gtable[i] = new TageEntry* [no_of_tiles];
	/* The above setup allows for a statically assigned max number of pointers to a tile. 
	   Now, changing a tile from one table to another is simply deallocating memory from one and allocating it in another. */
      }
	  DPRINTF(Tage, "\nBuilding TAGE Table of Pointers");
    }
    
    void
    DynamicTAGE::dynamicallyConfigTAGETable(std::vector<int> configuration_vector){
		for(int i = 1; i <= nHistoryTables; i++) {
			DPRINTF(Tage, "\n\n In History Table %d", i);
		//current_config is the config that was gotten after the previous dynamic TAGE config
		//thus current_config is the previous_config_vector
			uint8_t number_of_existing_tiles = current_config[i];
			DPRINTF(Tage, "\nNumber of Existing Tiles - %d", number_of_existing_tiles);
			DPRINTF(Tage, "\nConfig Vector Number of Tiles - %d", configuration_vector[i]);
		
			if(number_of_existing_tiles > configuration_vector[i]){
			//deallocate memory from those extra tiles
				for(int j = number_of_existing_tiles; j > configuration_vector[i]; j--){
					delete [] gtable[i][j];
					DPRINTF(Tage, "\n Deleting Tiles");
				}
			}
			
			else{
				for(int j = number_of_existing_tiles; j < configuration_vector[i]; j++){
					//Allocatre memory for new tiles
					gtable[i][j] = new TageEntry[entries_per_tile];
					DPRINTF(Tage, "\n Allocating Tiles to table");
				}
			}
		}
		current_config = configuration_vector;
		logTagTableSizes = current_config; //Do better with parameter passing
		DPRINTF(Tage, "\nMemory Alloc for TAGE Tables Complete");
	}
    
    
    bool DynamicTAGE::tagePredict(ThreadID tid, Addr branch_pc, bool cond_branch, BranchInfo* bi) {
	  DPRINTF(Tage, "\nIn TagePredict function");
      Addr pc = branch_pc;
      bool pred_taken = true;
      
      if (cond_branch){
		// TAGE prediction
		
		calculateIndicesAndTags(tid, pc, bi);
		
		bi->bimodalIndex = bindex(pc);
		
		bi->hitBank = 0;
		bi->hitTile = 0;
		bi->altBank = 0;
		
		int currentTableIndex = 0;
		int tileIndex = 0;
		int indexWithinTile = 0;
		
		//Look for the bank with longest matching history
		for (int i = nHistoryTables; i > 0; i--) {
			currentTableIndex = bi->tableIndices[i];
			tileIndex = calcTileIndex(currentTableIndex);
			indexWithinTile = calcIndexWithinTile(currentTableIndex);
			DPRINTF(Tage, "\nLook for the bank with longest matching history");
			if (noSkip[i] && gtable[i][tileIndex][indexWithinTile].tag == tableTags[i]) {
				bi->hitBank = i;
				bi->hitTile = tileIndex;
				bi->hitBankIndex = tableIndices[bi->hitBank];
				break;
			}
		
		}
		//Look for the alternate bank
		for (int i = bi->hitBank - 1; i > 0; i--) {
			DPRINTF(Tage, "\nLook for the alternate bank with longest matching history");
			currentTableIndex = bi->tableIndices[i];
			tileIndex = calcTileIndex(currentTableIndex);
			indexWithinTile = calcIndexWithinTile(currentTableIndex);
			
			if (noSkip[i] && gtable[i][tileIndex][indexWithinTile].tag == tableTags[i]) {
				bi->altBank = i;
				bi->altTile = tileIndex;
				bi->altBankIndex = tableIndices[bi->altBank];
				break;
			}
		}
		//computes the prediction and the alternate prediction
		if (bi->hitBank > 0) {
			if (bi->altBank > 0) {
				currentTableIndex = tableIndices[bi->altBank];
				tileIndex = calcTileIndex(currentTableIndex);
				indexWithinTile = calcIndexWithinTile(currentTableIndex);
				
				bi->altTaken =
				gtable[bi->altBank][tileIndex][indexWithinTile].ctr >= 0;
				extraAltCalc(bi);
			}
			else {
				bi->altTaken = getBimodePred(pc, bi);
			}
		
			currentTableIndex = tableIndices[bi->hitBank];
			tileIndex = calcTileIndex(currentTableIndex);
			indexWithinTile = calcIndexWithinTile(currentTableIndex);
			
			bi->longestMatchPred =
				gtable[bi->hitBank][tileIndex][indexWithinTile].ctr >= 0;
			
			currentTableIndex = bi->hitBankIndex;
			tileIndex = calcTileIndex(currentTableIndex);
			indexWithinTile = calcIndexWithinTile(currentTableIndex);
			
			bi->pseudoNewAlloc =
				abs(2 * gtable[bi->hitBank][tileIndex][indexWithinTile].ctr + 1) <= 1;
			
			//if the entry is recognized as a newly allocated entry and
			//useAltPredForNewlyAllocated is positive use the alternate
			//prediction
			if ((useAltPredForNewlyAllocated[getUseAltIdx(bi, branch_pc)] < 0)
				|| ! bi->pseudoNewAlloc) {
				bi->tagePred = bi->longestMatchPred;
				bi->provider = TAGE_LONGEST_MATCH;
			} else {
				bi->tagePred = bi->altTaken;
				bi->provider = bi->altBank ? TAGE_ALT_MATCH
				: BIMODAL_ALT_MATCH;
			}
		} 
		else {
			bi->altTaken = getBimodePred(pc, bi);
			bi->tagePred = bi->altTaken;
			bi->longestMatchPred = bi->altTaken;
			bi->provider = BIMODAL_ONLY;
		}
			//end TAGE prediction
		
		pred_taken = (bi->tagePred);
		DPRINTF(Tage, "Predict for %lx: taken?:%d, tagePred:%d, altPred:%d\n",
				branch_pc, pred_taken, bi->tagePred, bi->altTaken);
      }
      bi->branchPC = branch_pc;
      bi->condBranch = cond_branch;
      return pred_taken;
	}
	
    
    void DynamicTAGE::handleAllocAndUReset(bool alloc, bool taken, BranchInfo* bi, int nrand) {
		DPRINTF(Tage, "\nIn handleAllocAndUReset function");
      	if (alloc) {
			// is there some "unuseful" entry to allocate
			uint8_t min = 1;
			int currentTableIndex = 0;
			int tileIndex = 0;
			int indexWithinTile = 0;
			
			for (int i = nHistoryTables; i > bi->hitBank; i--) {
				currentTableIndex = bi->tableIndices[i];
				tileIndex = calcTileIndex(currentTableIndex);;
				indexWithinTile = calcIndexWithinTile(currentTableIndex);
				
				if (gtable[i][tileIndex][indexWithinTile].u < min) {
					min = gtable[i][tileIndex][indexWithinTile].u;
				}
			}
      
			// we allocate an entry with a longer history
			// to  avoid ping-pong, we do not choose systematically the next
			// entry, but among the 3 next entries
			int Y = nrand &
			((1ULL << (nHistoryTables - bi->hitBank - 1)) - 1);
			int X = bi->hitBank + 1;
			
			if (Y & 1) {
				X++;
				if (Y & 2)
					X++;
			}
      
			currentTableIndex = bi->tableIndices[X];
			tileIndex = calcTileIndex(currentTableIndex);;
			indexWithinTile = calcIndexWithinTile(currentTableIndex);	
			
			// No entry available, forces one to be available
			if (min > 0) {
				gtable[X][tileIndex][indexWithinTile].u = 0;
			}
      
      
			//Allocate entries
			unsigned numAllocated = 0;
			for (int i = X; i <= nHistoryTables; i++) {
				currentTableIndex = bi->tableIndices[i];
				tileIndex = calcTileIndex(currentTableIndex);;
				indexWithinTile = calcIndexWithinTile(currentTableIndex);	
				
				if (gtable[i][tileIndex][indexWithinTile].u == 0) {
					gtable[i][tileIndex][indexWithinTile].tag = bi->tableTags[i];
					gtable[i][tileIndex][indexWithinTile].ctr = (taken) ? 0 : -1;
					++numAllocated;
					if (numAllocated == maxNumAlloc) {
						break;
					}
				}
			}
      }
    
      tCounter++;
    
      handleUReset();
    }
  
    void DynamicTAGE::handleUReset() {
		DPRINTF(Tage, "\nIn handleUReset function");
      //periodic reset of u: reset is not complete but bit by bit
		if ((tCounter & ((1ULL << logUResetPeriod) - 1)) == 0) {
		// reset least significant bit
		// most significant bit becomes least significant bit
			for (int i = 1; i <= nHistoryTables; i++) {
				//for (int j = 0; j < (1ULL << logTagTableSizes[i]); j++) {
				for (int j = 0; j < current_config[i]; j++) {
					//accessing each tile
					for (int tileEntry = 0; tileEntry < entries_per_tile; tileEntry++) {
					//added this loop to access the tile entry		
					resetUctr(gtable[i][j][tileEntry].u);
					}
				}
			}
    	}
    }
  
  
    void DynamicTAGE::handleTAGEUpdate(Addr branch_pc, bool taken, BranchInfo* bi) {
		DPRINTF(Tage, "\nIn handleTAGEUpdate function");
  
		int currentTableIndex = 0;
		int tileIndex = 0;
		int indexWithinTile = 0;
	
		if (bi->hitBank > 0) {
			DPRINTF(Tage, "Updating tag table entry (%d,%d) for branch %lx\n",
			bi->hitBank, bi->hitBankIndex, branch_pc);
			
			currentTableIndex = bi->hitBankIndex;
			tileIndex = calcTileIndex(currentTableIndex);;
			indexWithinTile = calcIndexWithinTile(currentTableIndex);	
			
			ctrUpdate(gtable[bi->hitBank][tileIndex][indexWithinTile].ctr, taken, tagTableCounterBits);
			
			// if the provider entry is not certified to be useful also update
			// the alternate prediction
			if (gtable[bi->hitBank][tileIndex][indexWithinTile].u == 0) {
				if (bi->altBank > 0) {
					currentTableIndex = bi->altBankIndex;
					tileIndex = calcTileIndex(currentTableIndex);;
					indexWithinTile = calcIndexWithinTile(currentTableIndex);	
				
					ctrUpdate(gtable[bi->altBank][tileIndex][indexWithinTile].ctr, taken,
						tagTableCounterBits);
					DPRINTF(Tage, "Updating tag table entry (%d,%d) for"
						" branch %lx\n", bi->hitBank, bi->hitBankIndex,
						branch_pc);
				}
				if (bi->altBank == 0) {
					baseUpdate(branch_pc, taken, bi);
				}
			}
		
		
			currentTableIndex = bi->hitBankIndex;
			tileIndex = calcTileIndex(currentTableIndex);;
			indexWithinTile = calcIndexWithinTile(currentTableIndex);	
		
		
			// update the u counter
			if (bi->tagePred != bi->altTaken) {
			unsignedCtrUpdate(gtable[bi->hitBank][tileIndex][indexWithinTile].u, bi->tagePred == taken, tagTableUBits);
			}
    	} else {
		baseUpdate(branch_pc, taken, bi);
      }
    }



    // int8_t DynamicTAGE::getCtr(int hitBank, int hitBankIndex) const {
    //   int currentTableIndex = bi->hitBankIndex;
    //   tileIndex = calcTileIndex(currentTableIndex);;
    //   indexWithinTile = calcIndexWithinTile(currentTableIndex);	
  
    //   return gtable[hitBank][tileIndex][indexWithinTile].ctr;
    // }
    // //return speculativeHistUpdate;

  } // namespace branch_prediction
} // namespace gem5
