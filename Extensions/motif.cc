#include "motif.hh" 
#include "../out.hh" //Created when compiling Motif prediction algorithm.
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

static HashMap HairpinHashMap;
static HashMap InternalHashMap;
static HashMap BulgeHashMap;
static bool initialized;
static bool initialized2;
static bool initialized3;

std::vector<std::string> split(const std::string &s, char delim){
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;
    while (getline(ss,item,delim)) {
        result.push_back(item);
    }
    return result;   
}
HashMap globalMap_Hairpins(HashMap x) {
    std::string Line;
    std::ifstream infile("../RNAMotifs/Loops/HairpinMotifs.csv");
    while (getline (infile, Line)) {
        std::vector<std::string> v = split(Line,'+');
        std::string key = v[0];
        std::string prevalue = v[1];
        char value[prevalue.length()+1];
        strcpy(value,prevalue.c_str());
        char value2=value[0];
        x[key]=value2;
    }
    return x; 
}
HashMap globalMap_Internals(HashMap x){
    std::string Line;
    std::ifstream Data("../RNAMotifs/Loops/InternalMotifs.csv");
    while (getline (Data, Line)) {
        std::vector<std::string> v = split(Line,'+');
        std::string key = v[0];
        std::string prevalue = v[1];
        char value[prevalue.length()+1];
        strcpy(value,prevalue.c_str());
        char value2=value[0];
        x[key]=value2;
    }
    return x;
}
HashMap globalMap_Bulges(HashMap x){
    std::string Line;
    std::ifstream Data("../RNAMotifs/Loops/BulgeMotifs.csv");
    while (getline (Data, Line)) {
        std::vector<std::string> v = split(Line,'+');
        std::string key = v[0];
        std::string prevalue = v[1];
        char value[prevalue.length()+1];
        strcpy(value,prevalue.c_str());
        char value2=value[0];
        x[key]=value2;
    }
    return x;
}
std::string InputManagement(const Basic_Subsequence<char,unsigned int> &a){
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
    //std::cout << Motif << std::endl;
    return Motif; //Return the Motif string, correctly formatted for the HashMap
}
std::string InputManagement2(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b){
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
    //std::cout << Motif << std::endl;
    return Motif; //Return the Motif string, correctly formatted for the HashMap
}
char identify_motif(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '.';
    if (!initialized){
        initialized = true;
        HairpinHashMap = globalMap_Hairpins(HairpinHashMap);
    }
    std::string Motif;
    Motif = InputManagement(a); 
    if(HairpinHashMap.find(Motif) == HairpinHashMap.end()) {
        return res;
    }   else {
            res = HairpinHashMap[Motif];
            return res;
    }
}
char identify_motif(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b){
    char res = '.';
    if (!initialized2){
        initialized2 = true;
        InternalHashMap = globalMap_Internals(InternalHashMap);
    }   
    std::string Motif;
    Motif = InputManagement2(a,b); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    if(InternalHashMap.find(Motif) == InternalHashMap.end()) {
        return res;
    }
    else {
        res = InternalHashMap[Motif];
        return res;
    }
}
char identify_motif_b(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '.';
    if (!initialized3){
        initialized3 = true;
        BulgeHashMap = globalMap_Bulges(BulgeHashMap);
    }
    std::string Motif;
    Motif = InputManagement(a); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    if(BulgeHashMap.find(Motif) == BulgeHashMap.end()) {
        return res;
    }   else {
            res = BulgeHashMap[Motif];
            return res;
    }
}