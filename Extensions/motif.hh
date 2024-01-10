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
typedef std::unordered_map<std::string, char> HashMap;
static HashMap HairpinHashMap;
static HashMap InternalHashMap;
static HashMap BulgeHashMap;
static bool initialized;
static bool initialized2;
static bool initialized3;
static int Databa;
static int Databa2;
static int Databa3;
//This reads in the -Q arguemnt, choosing which database to use for motif predictions; 1 = BGSU, 2 = RMFam, 3 = both
inline static int Database() {
    return gapc::Opts::getOpts()->motifs;
}
//Split function that allows me to split my input strings from they xx+yy form into v[0]="xx" and v[1]="yy" splitting the sequences from the 
inline std::vector<std::string> split (const std::string &s, char delim){
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;
    while (getline(ss,item,delim)) {
        result.push_back(item);
    }
    return result;   
}
//HashMap implementation functions for all three loop types in the macrostate grammar.
inline HashMap globalMap_Hairpins(HashMap x, int DB) {
    std::string Line;
    const char* dab;
    if (DB == 1){dab = "Misc/Applications/RNAMotifs/Loops/HairpinMotifs.csv";}
        else if(DB == 2){dab ="Misc/Applications/RNAMotifs/Loops/HairpinRMFAM.csv";}
            else{dab ="Misc/Applications/RNAMotifs/Loops/HairpinALL.csv";}
    std::ifstream infile(dab);
    while (getline (infile, Line)) {
        std::vector<std::string> v = split(Line,'+');
        std::string key = v[0];
        char value[v[1].length()+1];
        strcpy(value,v[1].c_str());
        x[key]=value[0];
    }
    return x;
}
inline HashMap globalMap_Internals(HashMap x, int DB) {
    std::string Line;
    const char* dab2;
    if (DB == 1){dab2 = "Misc/Applications/RNAMotifs/Loops/InternalMotifs.csv";}
        else if(DB == 2){dab2 ="Misc/Applications/RNAMotifs/Loops/InternalRMFAM.csv";}
            else{dab2 ="Misc/Applications/RNAMotifs/Loops/InternalALL.csv";}
    std::ifstream Data2(dab2);
    while (getline (Data2, Line)) {
        std::vector<std::string> v = split(Line,'+');
        std::string key = v[0];
        char value[v[1].length()+1];
        strcpy(value,v[1].c_str());
        x[key]=value[0];
    }
    return x;
}
inline HashMap globalMap_Bulges(HashMap x, int DB) {
    std::string Line;
    const char* dab3;
    if (DB == 1){dab3 = "Misc/Applications/RNAMotifs/Loops/BulgeMotifs.csv";}
        else if(DB == 2){dab3 ="Misc/Applications/RNAMotifs/Loops/BulgeRMFAM.csv";}
            else{dab3 ="Misc/Applications/RNAMotifs/Loops/BulgeALL.csv";}
    std::ifstream Data3(dab3);
    while (getline (Data3, Line)) {
        std::vector<std::string> v = split(Line,'+');
        std::string key = v[0];
        char value[v[1].length()+1];
        strcpy(value,v[1].c_str());
        x[key]=value[0];
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
inline std::string InputManagement2(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b) {
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
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '.';
    if (!initialized){
        initialized = true;
        int Databa = Database();
        HairpinHashMap = globalMap_Hairpins(HairpinHashMap,Databa);
    }
    std::string Motif;
    Motif = InputManagement(a); 
    if(HairpinHashMap.find(Motif) == HairpinHashMap.end()) {
        return res;
    }   
    else {
            res = HairpinHashMap[Motif];
            return res;
    }
}
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b) {
    char res = '.';
    if (!initialized2){
        initialized2 = true;
        int Databa2 = Database();
        InternalHashMap = globalMap_Internals(InternalHashMap,Databa2);
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
inline char identify_motif_b(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '.';
    if (!initialized3){
        initialized3 = true;
        int Databa3 = Database();
        BulgeHashMap = globalMap_Bulges(BulgeHashMap,Databa3);
    }
    std::string Motif;
    Motif = InputManagement(a); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    //std::cout << Motif << std::endl;
    if(BulgeHashMap.find(Motif) == BulgeHashMap.end()) {
        return res;
    }   else {
            res = BulgeHashMap[Motif];
            return res;
    }
}
//Overloaded Motif function in shape version, only real difference here being that the base res character is '_' instead of '.' as this is the depiction used by the shape implementation
//the three _shape versions of the identify_motif functions are only used with shape levels 1 and 2 though, before that the shapes are too abstract to contain the extensions character '_'.
inline char identify_motif_shape(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '_';
    if (!initialized){
        initialized = true;
        int Databa = Database();
        HairpinHashMap = globalMap_Hairpins(HairpinHashMap,Databa);
    }
    std::string Motif;
    Motif = InputManagement(a); 
    if(HairpinHashMap.find(Motif) == HairpinHashMap.end()) {
        return res;
    }   
    else {
            res = HairpinHashMap[Motif];
            return res;
    }
}
inline char identify_motif_shape(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b) {
    char res = '_';
    if (!initialized2){
        initialized2 = true;
        int Databa2 = Database();
        InternalHashMap = globalMap_Internals(InternalHashMap,Databa2);
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
inline char identify_motif_b_shape(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '_';
    if (!initialized3){
        initialized3 = true;
        int Databa3 = Database();
        BulgeHashMap = globalMap_Bulges(BulgeHashMap,Databa3);
    }
    std::string Motif;
    Motif = InputManagement(a); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    //std::cout << Motif << std::endl;
    if(BulgeHashMap.find(Motif) == BulgeHashMap.end()) {
        return res;
    }   else {
            res = BulgeHashMap[Motif];
            return res;
    }
}
inline char identify_motif_hishape(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '\0';
    if (!initialized){
        initialized = true;
        int Databa = Database();
        HairpinHashMap = globalMap_Hairpins(HairpinHashMap,Databa);
    }
    std::string Motif;
    Motif = InputManagement(a); 
    if(HairpinHashMap.find(Motif) == HairpinHashMap.end()) {
        return res;
    }   
    else {
            res = HairpinHashMap[Motif];
            return res;
    }
}
inline char identify_motif_hishape(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b) {
    char res = '\0';
    if (!initialized2){
        initialized2 = true;
        int Databa2 = Database();
        InternalHashMap = globalMap_Internals(InternalHashMap,Databa2);
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
inline char identify_motif_b_hishape(const Basic_Subsequence<char, unsigned int> &a) {
    char res = '\0';
    if (!initialized3){
        initialized3 = true;
        int Databa3 = Database();
        BulgeHashMap = globalMap_Bulges(BulgeHashMap,Databa3);
    }
    std::string Motif;
    Motif = InputManagement(a); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    //std::cout << Motif << std::endl;
    if(BulgeHashMap.find(Motif) == BulgeHashMap.end()) {
        return res;
    }   else {
            res = BulgeHashMap[Motif];
            return res;
    }
}
#endif