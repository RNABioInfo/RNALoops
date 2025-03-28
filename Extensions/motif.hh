#pragma once
#include "subsequence.hh"
#include "rnaoptions.hh"
#include "mot_header.hh"
#include "shapes.hh"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <typeinfo>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <mutex>
#include <optional>
#include <tuple>
#include <algorithm>

struct Motif {std::string seq; char abb {'X'};};
struct directions {bool forward {false}; bool reverse {false};};
typedef std::unordered_map<std::string, char> HashMap;
static HashMap HairpinHashMap;
static HashMap InternalHashMap;
static HashMap BulgeHashMap;
static bool initialized;
static std::array Hairpins         = {rna3d_hairpins, rfam_hairpins, both_hairpins};
static std::array Hairpin_lengths  = {rna3d_hairpins_len, rfam_hairpins_len, both_hairpins_len};
static std::array Internals        = {rna3d_internals, rfam_internals, both_internals};
static std::array Internal_lengths = {rna3d_internals_len, rfam_internals_len, both_internals_len};
static std::array Bulges           = {rna3d_bulges, rfam_bulges, both_bulges};
static std::array Bulge_lengths    = {rna3d_bulges_len, rfam_bulges_len, both_bulges_len};

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

inline std::optional<Motif> parse_motif (const std::string &input, bool rev) {
    constexpr char delim = ',';
    std::stringstream sinput (input);
    size_t elements {0};
    Motif motif{};
    std::string item;
    while (std::getline(sinput,item,delim)) {
        if (elements == 0){
            if (rev){
                std::reverse(item.begin(),item.end());
            }
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
        std::cerr << "Found " << elements << ", Expected " << 2 << std::endl;
        return std::nullopt;
    }
    return motif;
}

inline std::string add_entry(HashMap &CurrentHashMap, std::string &line, bool reverse, std::string duplicate_string){
    duplicate_string = "";
    if (const auto motif = parse_motif(line, reverse)) {
       if (auto search = CurrentHashMap.find(motif.value().seq); search == CurrentHashMap.end()){
        CurrentHashMap[motif.value().seq]=motif.value().abb;
       }
       else{
           if (std::tolower(search->second) == std::tolower(motif.value().abb)){;}
           else{
            duplicate_string = motif.value().seq + "=" + motif.value().abb + "/" + CurrentHashMap[motif.value().seq];
            CurrentHashMap[motif.value().seq]=std::tolower(CurrentHashMap[motif.value().seq]);
            
           }
       }
    }
    return duplicate_string;
}

//HashMap implementation functions for all three loop types in the macrostate grammar.DB =database, which one gets used 1=BGSU, 2=RMFAM, 3 = bot // R = reverse, 1= No Reverses, 2 = Only Reverses, 3 = Both reverse and Forward
inline std::vector<std::string> Motif_HashMap(HashMap &CurrentHashMap, std::array<char* ,3> arr, std::array<unsigned int,3> len_arr, directions directions, std::vector<std::string> dupe_collector) {
    std::string str_mot_ar;
    switch(gapc::Opts::getOpts()->motifs) {
        case 1:
            str_mot_ar.assign(arr[0],len_arr[0]);
            break;
        case 2:
            str_mot_ar.assign(arr[1],len_arr[1]);
            break;
        case 3:
            str_mot_ar.assign(arr[2],len_arr[2]);
            break;
    }
    std::istringstream isstream(str_mot_ar);
    std::string line;
    std::string dupe;
    if (directions.forward) {
        while (std::getline (isstream, line,'\n')) {
            dupe = add_entry(CurrentHashMap,line,false,dupe);
            if (!dupe.empty()){
                dupe_collector.push_back(dupe);
            }
        }
    }
    isstream.clear();
    isstream.seekg(0);
    if (directions.reverse) {
        while (std::getline (isstream, line,'\n')) {
            dupe = add_entry(CurrentHashMap,line,true,dupe);
            if (!dupe.empty()) {
                dupe_collector.push_back(dupe);
            }
        }
    }
    return dupe_collector;
}

//Custom HashMap implementation function that fills HashMap X with the sequences saved in the csv file that path y points to. Capable of handling order reversing in 
 inline std::vector<std::string> Custom_Motif_HashMap(HashMap &CurrentHashMap, const std::string &path_to_csv, directions direction, std::vector<std::string> dupe_collector){
    std::string line;
    std::ifstream ifstream(path_to_csv);
    if (!ifstream.is_open()) {
        throw std::runtime_error("File not found");
    }
    std::string dupe;
    if (direction.forward){
        for (std::string line; std::getline(ifstream, line);) {
            dupe = add_entry(CurrentHashMap,line,false,dupe);
            if (!dupe.empty()) {
                dupe_collector.push_back(dupe);
            }
        }
    }
    ifstream.clear();
    ifstream.seekg(0);
    if (direction.reverse){
        for (std::string line; std::getline(ifstream, line);) {
            dupe = add_entry(CurrentHashMap,line,true,dupe);
            if (!dupe.empty()) {
                dupe_collector.push_back(dupe);
            }
        }
    }
    return dupe_collector;
}
//Input Manipulation Function, allowing for ONE RNA Basic_Subsequence inputs to be converted to the HashMap Key Formatting
inline std::string InputManagement(const Basic_Subsequence<char,unsigned int> &a) {
    std::string Motif;//Motif is initialized to later be the carrier of the actual sequence which is returned and later used to find the Motif in the HashMap
    for(const auto pos: a) {
        if (base_t(pos) == G_BASE) {
            Motif += "G";
        }   
        else if (base_t(pos) == A_BASE) {
            Motif += "A";
        }   
        else if (base_t(pos) == C_BASE) {
            Motif += "C";
        }   
        else if (base_t(pos) == U_BASE) {
            Motif += "U";
        }
    }
    return Motif; //Return the Motif string, correctly formatted for the HashMap
}

//Input Manipulation allowing for TWO RNA Basic_Subsequence inputs to be converted to the HashMap Key formatting
inline std::string InputManagement(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b) {
    std::string Motif1, Motif2, Motif;//Same as the Standard Version, except both sequences are converted to String first to then be connected via "$" sign
    for(const auto pos1: a) {
        if (base_t(pos1) == G_BASE) {
            Motif1 += "G";
        }   
        else if (base_t(pos1) == A_BASE) {
            Motif1 += "A";
        }   
        else if (base_t(pos1) == C_BASE) {
            Motif1 += "C";
        }   
        else if (base_t(pos1) == U_BASE) {
            Motif1 += "U";
        }
        throw std::runtime_error("Base of invalid type");
    }
    for(const auto pos2: b) {
        if (base_t(pos2) == G_BASE) {
            Motif += "G";
        }   
        else if (base_t(pos2) == A_BASE) {
            Motif += "A";
        }   
        else if (base_t(pos2) == C_BASE) {
            Motif += "C";
        }   
        else if (base_t(pos2) == U_BASE) {
            Motif += "U";
        }
        else {
            throw std::runtime_error("Base of invalid type");
        }

    }
    Motif = Motif1 + "$" + Motif2;
    return Motif; //Return the Motif string, correctly formatted for the HashMap
}

//Creates a Bellman's GAP style string from a Basic_Subsequence, this is different from the standard strings and not compatible with them!
inline String bgap_string(const Basic_Subsequence<char, unsigned int> &a) {
    String Motif;
    for (const auto sum: a) {
        if (base_t(sum) == G_BASE) {
            Motif.append('G'); //No += operator is implemented so we have to use append.
        }
        else if (base_t(sum) == A_BASE) {
            Motif.append('A');
        }
        else if (base_t(sum) == C_BASE) {
            Motif.append('C');
        }
        else {
            Motif.append('U');
        }
    }
    return Motif; // return the bgap style string
}

inline directions get_directions(directions &direction){
    switch(gapc::Opts::getOpts()->reversed){
        case 1:
            direction.forward = true;
            direction.reverse = false;
            break;
        case 2:
            direction.forward = false;
            direction.reverse = true;
            break;
        case 3:
            direction.forward = true;
            direction.reverse = true;
            break;
    }
    return direction;
}

//function that alters the global HashMaps (I know this might be bad practice but it's easy and it works.)
inline std::vector<std::string> fill_hashmap(std::string custom_path, bool custom_replace, HashMap &map, directions directions, std::array<char* ,3> arr, std::array<unsigned int,3> len_arr, std::vector<std::string> dupe_collector){
    if (!custom_path.empty()){
        if (custom_replace){
            dupe_collector = Custom_Motif_HashMap(map,custom_path, directions,  dupe_collector);
        }
        else{
            dupe_collector = Motif_HashMap(map,arr,len_arr,directions, dupe_collector);
            dupe_collector = Custom_Motif_HashMap(map,custom_path,directions, dupe_collector);
        }
    }
    else{
        dupe_collector = Motif_HashMap(map,arr,len_arr,directions,dupe_collector);
    }
    return dupe_collector;
}

inline void create_hashmaps(){
    std::vector<std::string> collector;
    directions directions;
    directions = get_directions(directions);
    const std::string& custom_hairpin_path = gapc::Opts::getOpts()->custom_hairpins;
    const std::string& custom_internal_path = gapc::Opts::getOpts()->custom_internals;
    const std::string& custom_bulge_path = gapc::Opts::getOpts() -> custom_bulges;
    bool replace_hairpins = gapc::Opts::getOpts() -> replaceH;
    bool replace_internals = gapc::Opts::getOpts() -> replaceI;
    bool replace_bulges = gapc::Opts::getOpts() -> replaceB;
    collector = fill_hashmap(custom_hairpin_path,  replace_hairpins,  HairpinHashMap,  directions, Hairpins,  Hairpin_lengths,  collector);
    collector = fill_hashmap(custom_internal_path, replace_internals, InternalHashMap, directions, Internals, Internal_lengths, collector);
    collector = fill_hashmap(custom_bulge_path,    replace_bulges,    BulgeHashMap,    directions, Bulges,    Bulge_lengths,    collector);
    if (collector.size() > 0){
    std::string dupe_string = std::accumulate(std::next(collector.begin()), collector.end(),collector[0],[](std::string a, std::string b){ return a +", " + b;});
    std::cerr << "Attention, duplicate sequences in motif set: " << dupe_string << "\n";}
    else{;}
}

//Overloaded identify_motif functions, two identify_motif for Hairpins and Internal Loops respectively while identify_motif_b is for bulge loops
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &a, char res) {
    if (!initialized){
        initialized = true;
        create_hashmaps();
    }
    std::string Motif;
    Motif = InputManagement(a);
    if (auto search = HairpinHashMap.find(Motif); search != HairpinHashMap.end()){
        return search->second;
    }
    return res;
}

//Experimental function that I want to later use for syntactic filtering, currently not in use
inline bool test_motif(const Basic_Subsequence<char, unsigned int> &a) {
    if (!initialized){
        initialized = true;
        create_hashmaps();
    }
    std::string Motif;
    Motif = InputManagement(a);
    if (auto search = HairpinHashMap.find(Motif); search != HairpinHashMap.end()){
        return true;
    }
    return false;
}

inline char identify_motif(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b, char res) {
    if (!initialized){
        initialized = true;
        create_hashmaps();
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
        create_hashmaps();
    }
    std::string Motif;
    Motif = InputManagement(a); //Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    if (auto search = BulgeHashMap.find(Motif); search != BulgeHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_align(const Basic_Subsequence<char, unsigned int> &a, const Basic_Subsequence<char, unsigned int> &b, char res) {
    if (!initialized) {
        initialized = true;
        create_hashmaps();
    }
    std::string Motif1;
    std::string Motif2;
    Motif1 = InputManagement(a); // Jedes mal wenn die Funktion aufgerufen wird, wird erst res erstellt und dann Input Management gerufen um die Basic_Subsequence zu verarbeiten
    Motif2 = InputManagement(b);
    char found1;
    char found2;
    if (auto search1 = HairpinHashMap.find(Motif1); search1 != HairpinHashMap.end()) {
        found1 = search1->second;
    }
    else {
        return res;
    }
    if (auto search2 = HairpinHashMap.find(Motif2); search2 != HairpinHashMap.end()) {
        found2 = search2->second;
    }
    else {
        return res;
    }
    if (std::toupper(found1) == std::toupper(found2)) {
        return found1;
    }
    else {
        return res;
    }
}

inline shape_t bl_shapeX(char mot, unsigned int shapelevel, shape_t &x){
    switch (shapelevel){
        case 5: 
            [[fallthrough]];
        case 4:
            if (mot != underScore) {
                return shape_t(mot) + x;
            }
            else{
                return x;
            }
        case 3:
            if (mot != underScore){
                return shape_t(openParen) + shape_t(mot) + x + shape_t(closeParen);
            }
            else{
                return shape_t(openParen) + x + shape_t(closeParen);
            }
        case 2:
            [[fallthrough]];
        case 1:
            return shape_t(openParen) + shape_t(mot) + x + shape_t(closeParen);
        default:
            std::cerr << "Shape level is not set" << std::endl;
            break;
    }
    throw std::invalid_argument("Shape level is not set");
}

inline shape_t br_shapeX(char mot, unsigned int shapelevel, shape_t &x){
    switch (shapelevel){
        case 5:
            [[fallthrough]];
        case 4:
            if (mot != underScore) {
                return x + shape_t(mot);
            }
            else{
                return x;
            }
        case 3:
            if (mot != underScore){
                return shape_t(openParen) + x + shape_t(mot) + shape_t(closeParen);
            }
            else{
                return shape_t(openParen) + x + shape_t(closeParen);
            }
        case 2:
            [[fallthrough]];
        case 1:
            return shape_t(openParen) + x + shape_t(mot) + shape_t(closeParen);
        default:
            std::cerr << "Shape level is not set" << std::endl;
            break;
    }
    throw std::invalid_argument("Shape level is not set");
}

inline shape_t il_shapeX(char mot, unsigned int shapelevel, shape_t &x){
    switch (shapelevel){
        case 5:
            if (mot != underScore) {
                return shape_t(mot) + x + shape_t(mot);
            }
            else{
                return x;
            }
        case 4:
            [[fallthrough]];
        case 3:
            if (mot != underScore) {
                return shape_t(openParen) + shape_t(mot) + x + shape_t(mot) + shape_t(closeParen);
            }
            else {
                return shape_t(openParen) + x + shape_t(closeParen);
            }
        case 2:
            [[fallthrough]];
        case 1:
            return shape_t(openParen) + shape_t(mot) + x + shape_t(mot) +shape_t(closeParen);
        default:
            std::cerr << "Shape level is not set" << std::endl;
            break;
    }
    throw std::invalid_argument("Shape level is not set");
}