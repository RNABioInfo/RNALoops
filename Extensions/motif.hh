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
#include "motifs/mot_header.hh"
typedef std::unordered_map<std::string, char> HashMap;
static HashMap HairpinHashMap;
static HashMap InternalHashMap;
static HashMap BulgeHashMap;
static bool initialized;
static std::array Hairpins         = {hbgsu_fw, hrfam_fw, hboth_fw, hbgsu_rv, hrfam_rv, hboth_rv, hbgsu_both, hrfam_both, hboth_both};
static std::array Hairpin_lengths  = {hbgsu_fw_len, hrfam_fw_len, hboth_fw_len, hbgsu_rv_len, hrfam_rv_len, hboth_rv_len, hbgsu_both_len, hrfam_both_len, hboth_both_len};
static std::array Internals        = {ibgsu_fw, irfam_fw, iboth_fw, ibgsu_rv, irfam_rv, iboth_rv, ibgsu_both, irfam_both, iboth_both};
static std::array Internal_lengths = {ibgsu_fw_len, irfam_fw_len, iboth_fw_len, ibgsu_rv_len, irfam_rv_len, iboth_rv_len, ibgsu_both_len, irfam_both_len, iboth_both_len};
static std::array Bulges           = {bbgsu_fw, brfam_fw, bboth_fw, bbgsu_rv, brfam_rv, bboth_rv, bbgsu_both, brfam_both, bboth_both};
static std::array Bulge_lengths    = {bbgsu_fw_len, brfam_fw_len, bboth_fw_len, bbgsu_rv_len, brfam_rv_len, bboth_rv_len, bbgsu_both_len, brfam_both_len, bboth_both_len};

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

//HashMap implementation functions for all three loop types in the macrostate grammar.DB =database, which one gets used 1=BGSU, 2=RMFAM, 3 = bot // R = reverse, 1= No Reverses, 2 = Only Reverses, 3 = Both reverse and Forward
inline HashMap Motif_HashMap(HashMap x, std::array<unsigned char* ,9> &arr, std::array<unsigned int,9> len_arr) {
    std::string str_mot_ar;
    switch(gapc::Opts::getOpts()->reversed){
       case 1:     
           switch(gapc::Opts::getOpts()->motifs){
               case 1:
                   str_mot_ar = reinterpret_cast<const char*>(arr[0]), len_arr[0];
                   break;
               case 2:
                   str_mot_ar = reinterpret_cast<const char*>(arr[1]), len_arr[1];
                   break;
               case 3:
                   str_mot_ar = reinterpret_cast<const char*>(arr[2]), len_arr[2];
                   break;
                   }
       break;
       case 2:
           switch(gapc::Opts::getOpts()->motifs){
               case 1:
                   str_mot_ar = reinterpret_cast<const char*>(arr[3]), len_arr[3];
                   break;
               case 2:
                   str_mot_ar = reinterpret_cast<const char*>(arr[4]), len_arr[4];
                   break;
               case 3:
                   str_mot_ar = reinterpret_cast<const char*>(arr[5]), len_arr[5];
                   break;
                   }
       break;
       case 3:
           switch(gapc::Opts::getOpts()->motifs){
                case 1:
                    str_mot_ar = reinterpret_cast<const char*>(arr[6]), len_arr[6];
                    break;
                case 2:
                    str_mot_ar = reinterpret_cast<const char*>(arr[7]), len_arr[7];
                    break;
                case 3:
                    str_mot_ar = reinterpret_cast<const char*>(arr[8]), len_arr[8];
                    break;
                   }
       break;
   }
   std::string Line;
   std::istringstream isstream(str_mot_ar);
   while (std::getline (isstream, Line,'\n')) {
       std::vector<std::string> v = split(Line,',');
       char value[v[1].length()];
       strcpy(value,v[1].c_str());
       x[v[0]]=value[0];
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
    else {
        return res;
    }
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
    else {
        return res;
    }
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
    else {
        return res;
    }
}
#endif