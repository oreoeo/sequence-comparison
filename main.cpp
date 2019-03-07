/* Bioinformatics Project 2
 *
 *
 *
 *
 */


#include <string>
#include <iostream>


// Constants
const int GAP_PENALTY = -2;
const int MISMATCH_PENALTY = -1;
const int MATCH_PENALTY = 1;

int main()
{
    // Initial strings of Shanghai and Ohio virus
    std::string seqX = "ECATCACCT";
    std::string seqY = "EGATACCC";
    
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
            std::cout << char_match;
            
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
            
//            std:: cout << matrix[i][j];
        }
        
//        std::cout << "\n";
    
    }
    
    // Find max score path in matrix
    
    
    std::cout << "\n";
    for(int i = 0; i < seqX.size(); i++)
    {
        for(int j = 0; j < seqY.size(); j++)
        {
            std:: cout << matrix[i][j];
        }
        std::cout << "\n";
    }
    
    
}
