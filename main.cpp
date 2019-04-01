/* Bioinformatics Project 2
 *
 *
 *
 *
 */


#include <string>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>


// Constants
const int GAP_PENALTY = -2;
const int MISMATCH_PENALTY = -3;
const int MATCH_PENALTY = 3;

// FIXME: How to keep from overflowing the stack when declaring our matrix in the program?
static int matrix[1507][1485];

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cerr << "You must supply filenames for the sequences when running the program. Terminating.\n";
    return 0;
  }
  else {
    std::ifstream seqXFile(argv[1]);
    std::ifstream seqYFile(argv[2]);

    // Initial strings of Shanghai and Ohio virus
    std::string seqX = "E";
    std::string seqY = "E";

    while (!seqXFile.eof()) {
      std::string tmp;
      seqXFile >> tmp;
      seqX += tmp;
    }

    while (!seqYFile.eof()) {
      std::string tmp;
      seqYFile >> tmp;
      seqY += tmp;
    }

    const int xSize = seqX.size();
    const int ySize = seqY.size();

    // Strings to hold the aligned sequences
    std::string alignedSeqX = "";
    std::string alignedSeqY = "";

    // Initial variables
    int max_score = 0;
    bool char_match = false;

    // Form a 2D array of our sequences
    // static int matrix[xSize][ySize];

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
    // TODO: Store locations of multiple absolute max scores for more than one possible solutions
    int absolute_max = matrix[0][0];
    int absolute_col = -1;
    int absolute_row = -1;

    // Navigate through the matrix finding the max values
    for (int j = 0; j < seqY.size(); ++j) {
      for (int i=0; i < seqX.size(); ++i) {
	// std::cout << matrix[i][j] << ' '; // Printing the matrix for testing
	if (matrix[i][j] > absolute_max) {
	  absolute_max = matrix[i][j];
	  absolute_col = i;
	  absolute_row = j;
	}
      }

      // std::cout << '\n'; // Printing the matrix for testing
    }

    int i = absolute_col;
    int j = absolute_row;
    
    // FIXME: Add gaps into the aligned sequences
    while (absolute_max > 0) {
      // DIAGONAL
      // Matching base
      if (absolute_max == matrix[i-1][j-1] + MATCH_PENALTY) {
	alignedSeqX += seqX[i];
	alignedSeqY += seqY[j];

	absolute_max = matrix[i-1][j-1];

	i -= 1;
	j -= 1;
      }
      // Mismatching base
      else if (absolute_max == matrix[i-1][j-1] + MISMATCH_PENALTY) {
	alignedSeqX += seqX[i];
	alignedSeqY += seqY[j];

	absolute_max = matrix[i-1][j-1];
	
	i -= 1;
	j -= 1;
      }

      // TOP GAP
      else if (absolute_max == matrix[i][j-1] + GAP_PENALTY) {
	break;
	alignedSeqX += "_";
	alignedSeqY += seqY[j];

	absolute_max = matrix[i][j-1];

	j -= 1;
      }

      // SIDE GAP
      else if (absolute_max == matrix[i-1][j] + GAP_PENALTY) {
	break;
	alignedSeqX += seqX[i];
	alignedSeqY += "_";

	absolute_max = matrix[i-1][j];
	
	i -= 1;
      }
    }

    reverse(alignedSeqX.begin(), alignedSeqX.end());
    reverse(alignedSeqY.begin(), alignedSeqY.end());


    std::cout << alignedSeqX << '\n' << alignedSeqY << '\n';
    std::cout << alignedSeqX.size() << ' ' << alignedSeqY.size() << '\n';

    return 0;
  }
}
