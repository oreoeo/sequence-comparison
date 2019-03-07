/* Bioinformatics Project 2
 *
 *
 *
 *
 */


#include <string>
#include <iostream>
#include <bits/stdc++.h>


// Constants
const int GAP_PENALTY = -2;
const int MISMATCH_PENALTY = -1;
const int MATCH_PENALTY = 1;

int main()
{
    // Initial strings of Shanghai and Ohio virus
    std::string seqX = "ECATCACCT";
    std::string seqY = "EGATACCC";

    std::string alignedSeqX = "";
    std::string alignedSeqY = "";
    
    // Initial variables
    int max_score = 0;
    bool char_match = false;
    
    // Form a 2D array of our sequences
    int matrix[seqX.size()][seqY.size()];
    
    // Init empty set
    for(int i = 0; i < seqX.size(); i++)
    {
        matrix[i][0] = 0;
    }
    for(int i = 0; i < seqY.size(); i++)
    {
        matrix[0][i] = 0;
    }
    
    // Populate matrix
    for(int i = 1; i < seqX.size(); i++)
    {
        for(int j = 1; j < seqY.size(); j++)
        {
            max_score = 0;
            char_match = false;
            
            if(seqX[i] == seqY[j])
            {
                char_match = true;
            }
            
            // Initially assume left case is the max score
            max_score = matrix[i-1][j] + GAP_PENALTY;
            
            // Check Diagonal score
            if(char_match)
            {
                if((matrix[i-1][j-1] + MATCH_PENALTY) > max_score)
                {
                    max_score = matrix[i-1][j-1] + MATCH_PENALTY;
                }
            }
            else
            {
                if ((matrix[i-1][j-1] + MISMATCH_PENALTY) > max_score)
                {
                    max_score = matrix[i-1][j-1] + MISMATCH_PENALTY;
                }
            }
            
            // Check top case
            if(matrix[i][j-1] + GAP_PENALTY > max_score)
            {
                max_score = matrix[i][j-1] + GAP_PENALTY;
            }
            
            // Check max score against 0
            if(max_score < 0)
            {
                max_score = 0;
            }
            
            // Set score in matrix
            matrix[i][j] = max_score;            
        }    
    }
    

    // Find the absolute max and where it is located
    int absolute_max = matrix[0][0];
    int absolute_col = -1;
    int absolute_row = -1;

    for(int i = 0; i < seqX.size(); i++)
    {
      for(int j = 0; j < seqY.size(); j++)
      {
    	if(matrix[i][j] > absolute_max)
    	{
    	  absolute_max = matrix[i][j];
    	  absolute_col = i;
    	  absolute_row = j;
	}
      }
    }

    // Trace back through the matrix and build our alignment
    int i = absolute_col;
    int j = absolute_row;
    
    while (true) {
	if (absolute_max == matrix[i-1][j] + GAP_PENALTY) {
    	  // Gap in left
    	  absolute_max -= matrix[i-1][j] + GAP_PENALTY;
	  i -= 1;
    	}
	
    	else if (absolute_max == matrix[i][j-1] + GAP_PENALTY) {
    	  // Gap in top
    	  absolute_max += GAP_PENALTY;
	  j -= 1;
    	}
	
    	else if ((absolute_max == matrix[i-1][j-1] + MATCH_PENALTY) || (absolute_max == matrix[i-1][j-1] + MISMATCH_PENALTY)) {
	  // Match the bases
	  alignedSeqX += seqX[i];
	  alignedSeqY += seqY[j];

    	  if (absolute_max == matrix[i-1][j-1] + MATCH_PENALTY) {
    	    absolute_max -= MATCH_PENALTY;
    	  }
    	  else if (absolute_max == matrix[i-1][j-1] + MISMATCH_PENALTY) {
    	    absolute_max += MISMATCH_PENALTY;
    	  }
	  i -= 1;
	  j -= 1;
    	}
	
	// Check if we are done
	if (absolute_max == 0) {
	  reverse(alignedSeqX.begin(), alignedSeqX.end());
	  reverse(alignedSeqY.begin(), alignedSeqY.end());
	  break;
	}
    }
    return 0;
}
