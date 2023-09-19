#ifndef MOTIF_HH
#define MOTIF_HH
#include "subsequence.hh"
typedef std::unordered_map<std::string, char> HashMap;
//Split function that allows me to split my input strings from they xx+yy form into v[0]="xx" and v[1]="yy" splitting the sequences from the 
std::vector<std::string> split (const std::string &s, char delim);
HashMap globalMap_Hairpins(HashMap x);
HashMap globalMap_Internals(HashMap x);
HashMap globalMap_Bulges(HashMap x);
//Input Manipulation Function, allowing for ONE RNA Basic_Subsequence inputs to be converted to the HashMap Key Formatting
std::string InputManagement(const Basic_Subsequence<char,unsigned int> &a);
//Input Manipulation allowing for TWO RNA Basic_Subsequence inputs to be converted to the HashMap Key formatting
std::string InputManagement2(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b);
//Overloaded identify_motif function, when one sequence comes in upper is used, when two are put in the lower is used
char identify_motif(const Basic_Subsequence<char, unsigned int> &a);
char identify_motif_b(const Basic_Subsequence<char, unsigned int> &a);
char identify_motif(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b);
#endif