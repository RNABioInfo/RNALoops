#ifndef MOTIF_HH
#define MOTIF_HH
#include "subsequence.hh"
#include "rnaoptions.hh"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<typeinfo>
#include<fstream>
#include<unordered_map>
#include<string>
#include<vector>
#include<sstream>
#include<map>
#include<mutex>
#include<optional>
#include "mot_header.hh"
struct Motif {std::string seq; char abb {'X'};};
typedef std::unordered_map<std::string, char> HashMap;
static HashMap HairpinHashMap;
static HashMap InternalHashMap;
static HashMap BulgeHashMap;
static bool initialized;
static std::array Hairpins         = {rna3d_hairpins_fw, rfam_hairpins_fw, both_hairpins_fw, rna3d_hairpins_rv, rfam_hairpins_rv, both_hairpins_rv, rna3d_hairpins_both, rfam_hairpins_both, both_hairpins_both,custom_hairpins};
static std::array Hairpin_lengths  = {rna3d_hairpins_fw_len, rfam_hairpins_fw_len, both_hairpins_fw_len, rna3d_hairpins_rv_len, rfam_hairpins_rv_len, both_hairpins_rv_len, rna3d_hairpins_both_len, rfam_hairpins_both_len, both_hairpins_both_len,custom_hairpins_len};
static std::array Internals        = {rna3d_internals_fw, rfam_internals_fw, both_internals_fw, rna3d_internals_rv, rfam_internals_rv, both_internals_rv, rna3d_internals_both, rfam_internals_both, both_internals_both,custom_internals};
static std::array Internal_lengths = {rna3d_internals_fw_len, rfam_internals_fw_len, both_internals_fw_len, rna3d_internals_rv_len, rfam_internals_rv_len, both_internals_rv_len, rna3d_internals_both_len, rfam_internals_both_len, both_internals_both_len,custom_internals_len};
static std::array Bulges           = {rna3d_bulges_fw, rfam_bulges_fw, both_bulges_fw, rna3d_bulges_rv, rfam_bulges_rv, both_bulges_rv, rna3d_bulges_both, rfam_bulges_both, both_bulges_both,custom_bulges};
static std::array Bulge_lengths    = {rna3d_bulges_fw_len, rfam_bulges_fw_len, both_bulges_fw_len, rna3d_bulges_rv_len, rfam_bulges_rv_len, both_bulges_rv_len, rna3d_bulges_both_len, rfam_bulges_both_len, both_bulges_both_len,custom_bulges_len};

//Split function that allows me to split my input strings from they xx,yy form into v[0]="xx" and v[1]="yy" splitting the sequences from their single letter abbreviations
inline std::vector<std::string> split (const std::string &s, char delim){
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;
    while (getline(ss,item,delim)) {
        result.push_back(item);
    }
    return result;   
}


inline std::optional<Motif> parse_motif (const std::string &input){
    constexpr char delim = ',';
    std::stringstream sinput (input);
    size_t elements {0};
    Motif motif{};
    std::string item;
    while (getline(sinput,item,delim)) {
        if (elements == 0){
            motif.seq = item;
        }
        else if (elements == 1){
            motif.abb = item.front();
        }
        else{
            std::cerr << "Error more than one motif in one line" << std::endl;
            return std::nullopt;
        }
        ++elements;
    }
    if (elements > 2){
        std::cerr << "Found " << elements << ", Expecetd " << 2 << std::endl;
        return std::nullopt;
    }
    return motif;
}


//HashMap implementation functions for all three loop types in the macrostate grammar.DB =database, which one gets used 1=BGSU, 2=RMFAM, 3 = bot // R = reverse, 1= No Reverses, 2 = Only Reverses, 3 = Both reverse and Forward
inline HashMap Motif_HashMap(HashMap x, std::array<char* ,10> arr, std::array<unsigned int,10> len_arr) {
    std::string str_mot_ar;
    switch(gapc::Opts::getOpts()->reversed){
       case 1:     
           switch(gapc::Opts::getOpts()->motifs){
               case 1:
                   str_mot_ar.assign(arr[0],len_arr[0]);
                   break;
               case 2:
                   str_mot_ar.assign(arr[1],len_arr[1]);
                   break;
               case 3:
                   str_mot_ar.assign(arr[2],len_arr[2]);
                   break;
                case 4:
                   str_mot_ar.assign(arr[9],len_arr[9]);
                   break;
                   }
       break;
       case 2:
           switch(gapc::Opts::getOpts()->motifs){
               case 1:
                   str_mot_ar.assign(arr[3],len_arr[3]);
                   break;
               case 2:
                   str_mot_ar.assign(arr[4],len_arr[4]);
                   break;
               case 3:
                   str_mot_ar.assign(arr[5],len_arr[5]);
                   break;
                case 4:
                   str_mot_ar.assign(arr[9],len_arr[9]);
                   break;
                   }
       break;
       case 3:
           switch(gapc::Opts::getOpts()->motifs){
                case 1:
                    str_mot_ar.assign(arr[6],len_arr[6]);
                    break;
                case 2:
                    str_mot_ar.assign(arr[7],len_arr[7]);
                    break;
                case 3:
                    str_mot_ar.assign(arr[8],len_arr[8]);
                    break;
                case 4:
                    str_mot_ar.assign(arr[9],len_arr[9]);
                    break;
                   }
       break;
   }
   std::string line;
   std::istringstream isstream(str_mot_ar);
   while (std::getline (isstream, line,'\n')) {
        if (const auto motif = parse_motif(line)){
            x[motif.value().seq]=motif.value().abb;
        }
   }
    return x;
}

//Input Manipulation Function, allowing for ONE RNA Basic_Subsequence inputs to be converted to the HashMap Key Formatting
inline std::string InputManagement(const Basic_Subsequence<char,unsigned int> &a) {
    std::string Motif;//Motif is initialized to later be the carrier of the actual sequence which is returned and later used to find the Motif in the HashMap
    for(unsigned int p = 0; p < a.size();p++){
        int sum;
        sum = a.i + p;
        if (base_t(a.seq->seq[sum]) == G_BASE){
            Motif += "G";
        }   else if (base_t(a.seq->seq[sum]) == A_BASE){
            Motif += "A";
        }   else if (base_t(a.seq->seq[sum]) == C_BASE){
            Motif += "C";
        }   else {
            Motif += "U";
        }
    }
    return Motif; //Return the Motif string, correctly formatted for the HashMap
}

//Input Manipulation allowing for TWO RNA Basic_Subsequence inputs to be converted to the HashMap Key formatting
inline std::string InputManagement(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b) {
    std::string Motif1, Motif2, Motif;//Same as the Standard Version, except both sequences are converted to String first to then be connected via "$" sign
    for(unsigned int t = 0; t < a.size();t++){
        int sum;
        sum = a.i + t;
        if (base_t(a.seq->seq[sum]) == G_BASE){
            Motif1 += "G";           
        }   else if (base_t(a.seq->seq[sum]) == A_BASE){
            Motif1 += "A";
        }   else if (base_t(a.seq->seq[sum]) == C_BASE){
            Motif1 += "C";
        }   else {
            Motif1 += "U";
        }
    }
    for(unsigned int t2 = 0; t2 < b.size();t2++){
        int sum;
        sum = b.i + t2;
        if (base_t(b.seq->seq[sum]) == G_BASE){
            Motif2 += "G";
        }   else if (base_t(b.seq->seq[sum]) == A_BASE){
            Motif2 += "A";
        }   else if (base_t(b.seq->seq[sum]) == C_BASE){
            Motif2 += "C";
        }   else {
            Motif2 += "U";
        }
    }
    Motif = Motif1 + "$" + Motif2;
    return Motif; //Return the Motif string, correctly formatted for the HashMap
}

//Overloaded identify_motif functions, two identify_motif for Hairpins and Internal Loops respectively while identify_motif_b is for bulge loops
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &a, char res) {
    if (!initialized){
        initialized = true;
        HairpinHashMap  = Motif_HashMap(HairpinHashMap, Hairpins, Hairpin_lengths);
        InternalHashMap = Motif_HashMap(InternalHashMap, Internals, Internal_lengths);
        BulgeHashMap    = Motif_HashMap(BulgeHashMap, Bulges, Bulge_lengths);
    }
    std::string Motif;
    Motif = InputManagement(a);
    if (auto search = HairpinHashMap.find(Motif); search != HairpinHashMap.end()){
        return search->second;
    }
    return res;
}

inline char identify_motif(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b, char res) {
    if (!initialized){
        initialized = true;
        HairpinHashMap  = Motif_HashMap(HairpinHashMap, Hairpins, Hairpin_lengths);
        InternalHashMap = Motif_HashMap(InternalHashMap, Internals, Internal_lengths);
        BulgeHashMap    = Motif_HashMap(BulgeHashMap, Bulges, Bulge_lengths);
    }   
    std::string Motif;
    Motif = InputManagement(a,b); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    if (auto search = InternalHashMap.find(Motif); search != InternalHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_b(const Basic_Subsequence<char, unsigned int> &a, char res) {
    if (!initialized){
        initialized = true;
        HairpinHashMap  = Motif_HashMap(HairpinHashMap, Hairpins, Hairpin_lengths);
        InternalHashMap = Motif_HashMap(InternalHashMap, Internals, Internal_lengths);
        BulgeHashMap    = Motif_HashMap(BulgeHashMap, Bulges, Bulge_lengths);
    }
    std::string Motif;
    Motif = InputManagement(a); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    if (auto search = BulgeHashMap.find(Motif); search != BulgeHashMap.end()) {
        return search->second;
    }
    return res;
}
#endif