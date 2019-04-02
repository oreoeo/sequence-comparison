/* Bioinformatics Project 2
 *
 *
 *
 *
 */


#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <bits/stdc++.h>

std::string transcribeToAminoAcidSequence(const std::string&);
int countIndels(const std::string&, const std::string&);
int countSynMutations(const std::string&, const std::string&);
int countNonSynMutations(const std::string&, const std::string&);

// Constants
const int GAP_PENALTY = -2;
const int MISMATCH_PENALTY = -1;
const int MATCH_PENALTY = 1;

const std::map<std::string, std::string> AMINO_ACIDS {
						      {"UUU","Phe"},
						      {"UUC","Phe"},
						      {"UUA","Leu"},
						      {"UUG","Leu"},
						      {"CUU","Leu"},
						      {"CUC","Leu"},
						      {"CUA","Leu"},
						      {"CUG","Leu"},
						      {"UCU","Ser"},
						      {"UCC","Ser"},
						      {"UCA","Ser"},
						      {"UCG","Ser"},
						      {"UAU","Tyr"},
						      {"UAC","Tyr"},
						      {"UAA","Stp"},
						      {"UAG","Stp"},
						      {"UGU","Cys"},
						      {"UGC","Cys"},
						      {"UGA","Stp"},
						      {"UGG","Trp"},
						      {"CCU","Pro"},
						      {"CCC","Pro"},
						      {"CCA","Pro"},
						      {"CCG","Pro"},
						      {"CAU","His"},
						      {"CAC","His"},
						      {"CAA","Gln"},
						      {"CAG","Gln"},
						      {"CGU","Arg"},
						      {"CGC","Arg"},
						      {"CGA","Arg"},
						      {"CGG","Arg"},
						      {"AUU","Ile"},
						      {"AUC","Ile"},
						      {"AUA","Ile"},
						      {"AUG","Met"},
						      {"ACU","Thr"},
						      {"ACC","Thr"},
						      {"ACA","Thr"},
						      {"ACG","Thr"},
						      {"AAU","Asn"},
						      {"AAC","Asn"},
						      {"AAA","Lys"},
						      {"AAG","Lys"},
						      {"AGU","Ser"},
						      {"AGC","Ser"},
						      {"AGA","Arg"},
						      {"AGG","Arg"},
						      {"GUU","Val"},
						      {"GUC","Val"},
						      {"GUA","Val"},
						      {"GUG","Val"},
						      {"GCU","Ala"},
						      {"GCC","Ala"},
						      {"GCA","Ala"},
						      {"GCG","Ala"},
						      {"GAU","Asp"},
						      {"GAC","Asp"},
						      {"GAA","Glu"},
						      {"GAG","Glu"},
						      {"GGU","Gly"},
						      {"GGC","Gly"},
						      {"GGA","Gly"},
						      {"GGG","Gly"}
};

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
    std::ifstream blastQFile("blast_qseq.txt");
    std::ifstream blastHFile("blast_hseq.txt");

    // Initial strings of Shanghai and Ohio virus
    std::string seqX = "E";
    std::string seqY = "E";
    std::string blastQ = "";
    std::string blastH = "";

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

    while (!blastQFile.eof()) {
      std::string tmp;
      blastQFile >> tmp;
      blastQ += tmp;
    }

    while (!blastHFile.eof()) {
      std::string tmp;
      blastHFile >> tmp;
      blastH += tmp;
    }
    
    const int xSize = seqX.size();
    const int ySize = seqY.size();

    std::vector<std::pair<std::string,std::string>> aminoAcidSequences;

    // Front of the vector will always be the original amino acid sequences
    aminoAcidSequences.push_back(std::pair<std::string,std::string>(transcribeToAminoAcidSequence(seqX),transcribeToAminoAcidSequence(seqY)));
    aminoAcidSequences.push_back(std::pair<std::string,std::string>(transcribeToAminoAcidSequence(blastQ),transcribeToAminoAcidSequence(blastH)));    

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
    int absolute_max = matrix[0][0];
    int absolute_col = -1;
    int absolute_row = -1;

    std::vector<std::pair<int,int>> max_score_locations;

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

    // Store the absolute max score we found in our matrix
    // Check for other occurences of the max score in the matrix
    max_score_locations.push_back(std::pair<int, int>(absolute_col, absolute_row));

    for (int j = 0; j < seqY.size(); ++j) {
      for (int i = 0; i < seqX.size(); ++i) {
	if ((matrix[i][j] == absolute_max) &&
	    (max_score_locations.front().first != i) && (max_score_locations.front().first != i)) {
	  max_score_locations.push_back(std::pair<int, int>(i, j));
	}
      }
    }

    // Hold our solution(s) and trace back to find what they are
    std::vector<std::pair<std::string,std::string>> alignedSequences;

    for (int i = 0; i < max_score_locations.size(); ++i) {
      int col = max_score_locations[i].first;
      int row = max_score_locations[i].second;

      int score = absolute_max;
      alignedSeqX = "";
      alignedSeqY = "";

      while (score > 0) {
	// DIAGONAL
	// Matching base
	if (score == matrix[col-1][row-1] + MATCH_PENALTY) {
	  alignedSeqX += seqX[col];
	  alignedSeqY += seqY[row];

	  score = matrix[col-1][row-1];

	  col -= 1;
	  row -= 1;
	}
	// Mismatching base
	else if (score == matrix[col-1][row-1] + MISMATCH_PENALTY) {
	  alignedSeqX += seqX[col];
	  alignedSeqY += seqY[row];

	  score = matrix[col-1][row-1];
	
	  col -= 1;
	  row -= 1;
	}

	// TOP GAP
	else if (score == matrix[col][row-1] + GAP_PENALTY) {
	  alignedSeqX += "_";
	  alignedSeqY += seqY[row];

	  score = matrix[col][row-1];

	  row -= 1;
	}

	// SIDE GAP
	else if (score == matrix[col-1][row] + GAP_PENALTY) {
	  alignedSeqX += seqX[col];
	  alignedSeqY += "_";

	  score = matrix[col-1][row];
	
	  col -= 1;
	}
      }
      
      reverse(alignedSeqX.begin(), alignedSeqX.end());
      reverse(alignedSeqY.begin(), alignedSeqY.end());

      alignedSequences.push_back(std::pair<std::string,std::string>(alignedSeqX, alignedSeqY));

      aminoAcidSequences.push_back(std::pair<std::string,std::string>(transcribeToAminoAcidSequence(alignedSeqX),transcribeToAminoAcidSequence(alignedSeqY)));
    }

    // Output our alignment findings
    std::ofstream outFile("aligned_sequences.txt");

    for (int i = 0; i < alignedSequences.size(); ++i) {
      outFile << "Alignment " << i << '\n';
      outFile << "SeqX: " << alignedSequences[i].first << '\n';
      outFile << "SeqY: " << alignedSequences[i].second << '\n';
      outFile << "Length of alignment: " << alignedSequences[i].first.size() << '\n';
      outFile << "--------------------\n";

      outFile << "Original SeqX amino acid sequence: " << aminoAcidSequences[0].first << '\n';
      outFile << "Aligned SeqX amino acid sequence: " << aminoAcidSequences[i+2].first << '\n';
      outFile << "SeqX indels: " << countIndels(aminoAcidSequences[0].first, aminoAcidSequences[i+2].first) << '\n';
      outFile << "SeqX synonymous mutations: " << countSynMutations(aminoAcidSequences[0].first, aminoAcidSequences[i+2].first) << '\n';
      outFile << "SeqX non-synonymous mutations: " << countNonSynMutations(aminoAcidSequences[0].first, aminoAcidSequences[i+2].first) << '\n';
      outFile << "--------------------\n";


      outFile << "Original SeqY amino acid sequence: " << aminoAcidSequences[0].second << '\n';
      outFile << "Aligned SeqY amino acid sequence: " << aminoAcidSequences[i+2].second << '\n';
      outFile << "SeqY indels: " << countIndels(aminoAcidSequences[0].second, aminoAcidSequences[i+2].second) << '\n';
      outFile << "SeqY synonymous mutations: " << countSynMutations(aminoAcidSequences[0].second, aminoAcidSequences[i+2].second) << '\n';
      outFile << "SeqY non-synonymous mutations: " << countNonSynMutations(aminoAcidSequences[0].second, aminoAcidSequences[i+2].second) << '\n';
      outFile << "====================\n";
    }

    outFile << "BLAST Results\n";
    outFile << "BlastQ: " << blastQ << '\n';
    outFile << "BlastH: " << blastH << '\n';
    outFile << "--------------------\n";

    outFile << "BlastQ indels: " << countIndels(aminoAcidSequences[0].first, aminoAcidSequences[1].first) << '\n';
    outFile << "BlastQ synonymous mutations: " << countSynMutations(aminoAcidSequences[0].first, aminoAcidSequences[1].first) << '\n';
    outFile << "BlastQ non-synonymous mutations: " << countNonSynMutations(aminoAcidSequences[0].first, aminoAcidSequences[1].first) << '\n';
    outFile << "--------------------\n";

    outFile << "BlastH indels: " << countIndels(aminoAcidSequences[0].second, aminoAcidSequences[1].second) << '\n';
    outFile << "BlastH synonymous mutations: " << countSynMutations(aminoAcidSequences[0].second, aminoAcidSequences[1].second) << '\n';
    outFile << "BlastH non-synonymous mutations: " << countNonSynMutations(aminoAcidSequences[0].second, aminoAcidSequences[1].second) << '\n';
    outFile << "====================\n";

    return 0;
  }
}

std::string transcribeToAminoAcidSequence(const std::string& nuc_seq) {
  std::string aa_seq = "";

  for (int i = 0; i+3 < nuc_seq.size(); i+=3)
    aa_seq += AMINO_ACIDS.find(nuc_seq.substr(i,3))->second;

  return aa_seq;
}

int countIndels(const std::string& aaSeqX, const std::string& aaSeqY) {
  if (aaSeqX.size() >= aaSeqY.size())
    return aaSeqX.size() - aaSeqY.size();
  
  else if (aaSeqY.size() > aaSeqX.size()) 
    return aaSeqY.size() - aaSeqX.size();

  return -1;
}

int countSynMutations(const std::string& aaSeqX, const std::string& aaSeqY) {
  int synMutation = 0;

  int smallerSize;

  if (aaSeqX.size() <= aaSeqY.size())
    smallerSize = aaSeqX.size();
  else
    smallerSize = aaSeqY.size();

  for (int i = 0; i+3 < smallerSize; i+=3)
    if (aaSeqX.substr(i,3) == aaSeqY.substr(i,3))
      synMutation += 1;
  
  return synMutation;
}

int countNonSynMutations(const std::string& aaSeqX, const std::string& aaSeqY) {
  int nonSynMutation = 0;

  int smallerSize;

  if (aaSeqX.size() <= aaSeqY.size())
    smallerSize = aaSeqX.size();
  else
    smallerSize = aaSeqY.size();

  for (int i = 0; i+3 < smallerSize; i+=3)
    if (aaSeqX.substr(i,3) != aaSeqY.substr(i,3))
      nonSynMutation += 1;
  
  return nonSynMutation;
}
