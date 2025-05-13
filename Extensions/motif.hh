#pragma once
#include "ali_t.hh"
#include "rna.hh"
#include "rnaoptions_defaults.hh"
#include "rope.hh"
#include "subsequence.hh"
#include "sequence.hh"
#include "rnaoptions.hh"
#include "mot_header.hh"
#include "shapes.hh"
#include "shape.hh"
#include <cctype>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sstream>
#include <optional>
#include <algorithm>

struct size2_vector_hash {
    template <class T1>
    size_t operator () (const std::vector<T1> &input_vector) const {
        auto value1 = std::hash<T1>{}(input_vector[0]);
        auto value2 = std::hash<T1>{}(input_vector[1]);
        return value1 ^( value2 << 1);
    }
};

using shape_t =  Shape;
struct Motif {std::string seq; char abb {'X'};};
struct directions {bool forward {false}; bool reverse {false};};
using HashMap = std::unordered_map<Basic_Sequence<char, unsigned int>, char, Hash_ali_array>;
static HashMap HairpinHashMap;
static HashMap InternalHashMap;
static HashMap BulgeHashMap;
static std::unordered_map<std::vector<char>, char, size2_vector_hash> DupeHashMap;
static bool initialized;
static std::array Hairpins         = {rna3d_hairpins, rfam_hairpins};
static std::array Hairpin_lengths  = {rna3d_hairpins_len, rfam_hairpins_len};
static std::array Internals        = {rna3d_internals, rfam_internals};
static std::array Internal_lengths = {rna3d_internals_len, rfam_internals_len,};
static std::array Bulges           = {rna3d_bulges, rfam_bulges};
static std::array Bulge_lengths    = {rna3d_bulges_len, rfam_bulges_len};

enum shapelevel_enum: std::uint8_t {five=5,four=4,three=3,two=2,one=1};

//Split function that allows me to split my input strings from they xx,yy form into v[0]="xx" and v[1]="yy" splitting the sequences from their single letter abbreviations
inline std::vector<std::string> split (const std::string &input_str, char delim){
    std::vector<std::string> result;
    std::stringstream streamd_string (input_str);
    std::string item;
    while (std::getline(streamd_string,item,delim)) {
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
            std::cerr << "Error more than one motif in one line" << "\n";
            return std::nullopt;
        }
        ++elements;
    }
    if (elements > 2){
        std::cerr << "Found " << elements << ", Expected " << 2 << "\n";
        return std::nullopt;
    }
    return motif;
}

inline void add_entry(HashMap &CurrentHashMap, std::string &line, bool reverse, std::vector<char>& dupe_collector){
    if (const auto motif = parse_motif(line, reverse)) {
        Basic_Sequence<char, unsigned int> basic_seq_motif {motif.value().seq.data(), static_cast<unsigned int>(motif.value().seq.size())};
        char_to_rnali(basic_seq_motif);
       if (auto search = CurrentHashMap.find(basic_seq_motif); search == CurrentHashMap.end()){
        CurrentHashMap[basic_seq_motif]=motif.value().abb;
       }
       else{
           if (std::tolower(search->second) == std::tolower(motif.value().abb,std::locale())){;}
           else{//This is where we go if the same sequences has different motifs

                std::vector<char> dupes {std::toupper(search->second,std::locale()),std::toupper(motif.value().abb,std::locale())};

                std::sort(dupes.begin(),dupes.end());
                if (auto dupe_search = DupeHashMap.find(dupes); dupe_search == DupeHashMap.end()){ //If these two motifs haven't appeared as a duplicate yet:
                
                    if (std::find(dupe_collector.begin(),dupe_collector.end(),std::tolower(search->second)) == dupe_collector.end()){ // If the first element of dupes is not used yet:
                        CurrentHashMap[basic_seq_motif] = std::tolower(search->second,std::locale());
                        DupeHashMap[dupes] = std::tolower(search->second,std::locale());
                        dupe_collector.push_back(std::tolower(search->second,std::locale()));
                    }

                    else if (std::find(dupe_collector.begin(),dupe_collector.end(),std::tolower(motif.value().abb)) == dupe_collector.end()){ // If the second element of dupes is also not used yet:
                            CurrentHashMap[basic_seq_motif] = std::tolower(motif.value().abb,std::locale());
                            DupeHashMap[dupes] = std::tolower(motif.value().abb,std::locale());
                            dupe_collector.push_back(std::tolower(motif.value().abb,std::locale()));
                    }

                    else { //Threeway sequence, throw errro cause this aint implemented (yet)
                        throw std::runtime_error("A single sequence is bound to three different motif, this WILL lead to issues. No fix for this currently available.");
                    }
                }  
                else {
                    CurrentHashMap[basic_seq_motif] = dupe_search->second;
            }
           }
       }
    }
}

//HashMap implementation functions for all three loop types in the macrostate grammar.DB =database, which one gets used 1=BGSU, 2=RMFAM, 3 = both // R = reverse, 1= No Reverses, 2 = Only Reverses, 3 = Both reverse and Forward
inline std::vector<char> Motif_HashMap(HashMap &CurrentHashMap, std::array<unsigned char* ,2> arr, std::array<unsigned int,2> len_arr, directions directions, std::vector<char> dupe_collector) {
    std::string str_mot_ar;
    switch(gapc::Opts::getOpts()->motifs) {
        case 1:
            str_mot_ar.append(reinterpret_cast<char*>(arr[0]),len_arr[0]);
            break;
        case 2:
            str_mot_ar.append(reinterpret_cast<char*>(arr[1]),len_arr[1]);
            break;
        case 3:
            str_mot_ar.append(reinterpret_cast<char*>(arr[0]),len_arr[0]);
            str_mot_ar.append(reinterpret_cast<char*>(arr[1]),len_arr[1]);
            break;
        default:
            throw std::runtime_error("Could not identify given built in motif sets");
    }
    std::istringstream isstream(str_mot_ar);
    std::string line;
    if (directions.forward) {
        while (std::getline (isstream, line,'\n')) {
            add_entry(CurrentHashMap,line,false,dupe_collector);
        }
    }
    isstream.clear();
    isstream.seekg(0);
    if (directions.reverse) {
        while (std::getline (isstream, line,'\n')) {
            add_entry(CurrentHashMap,line,true,dupe_collector);
        }
    }
    return dupe_collector;
}

//Custom HashMap implementation function that fills HashMap X with the sequences saved in the csv file that path y points to. Capable of handling order reversing
 inline std::vector<char> Custom_Motif_HashMap(HashMap &CurrentHashMap, const std::string &path_to_csv, directions direction, std::vector<char> dupe_collector){
    std::ifstream ifstream(path_to_csv);
    if (!ifstream.is_open()) {
        throw std::runtime_error("File not found");
    }
    if (direction.forward){
        for (std::string line; std::getline(ifstream, line);) {
            add_entry(CurrentHashMap,line,false,dupe_collector);
        }
    }
    ifstream.clear();
    ifstream.seekg(0);
    if (direction.reverse){
        for (std::string line; std::getline(ifstream, line);) {
            add_entry(CurrentHashMap,line,true,dupe_collector);
        }
    }
    return dupe_collector;
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
        default:
            throw std::runtime_error("Something went wrong with the directions...");
    }
    return direction;
}

//function that alters the global HashMaps (I know this might be bad practice but it's easy and it works.)
inline void fill_hashmap(const std::string& custom_path, bool custom_replace, HashMap &map, directions directions, std::array<unsigned char* ,2> arr, std::array<unsigned int,2> len_arr, std::vector<char> &dupe_collector){
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
}
// 2025-05-05 16:26:27:results:WARNING: Attention, duplicate sequences in motif set: UGAGAAU=G/T -> g, GGAA=G/U -> g, GUGA=U/G -> u

inline void create_hashmaps(){
    std::vector<char> collector;
    directions directions;
    directions = get_directions(directions);
    const std::string& custom_hairpin_path = gapc::Opts::getOpts()->custom_hairpins;
    const std::string& custom_internal_path = gapc::Opts::getOpts()->custom_internals;
    const std::string& custom_bulge_path = gapc::Opts::getOpts() -> custom_bulges;
    bool replace_hairpins = gapc::Opts::getOpts() -> replaceH;
    bool replace_internals = gapc::Opts::getOpts() -> replaceI;
    bool replace_bulges = gapc::Opts::getOpts() -> replaceB;
    fill_hashmap(custom_hairpin_path,  replace_hairpins,  HairpinHashMap,  directions, Hairpins,  Hairpin_lengths,  collector);
    fill_hashmap(custom_internal_path, replace_internals, InternalHashMap, directions, Internals, Internal_lengths, collector);
    fill_hashmap(custom_bulge_path,    replace_bulges,    BulgeHashMap,    directions, Bulges,    Bulge_lengths,    collector);
    std::vector<std::string> strings;
    for (const auto &keyvalue: DupeHashMap){
        std::string this_kv = std::string() + keyvalue.first[0] + "/" + keyvalue.first[1] + " -> " + keyvalue.second;
        strings.emplace_back(this_kv);
        }
        std::string dupe_string = std::accumulate(std::next(strings.begin()), strings.end(),strings[0],[](const std::string& Existing_string, const std::string& Added_string){ return Existing_string +", " + Added_string;});
        std::cerr << "Attention, the same sequences appears for different motifs. Ambiguity cases: " << dupe_string << "\n";
}

//Overloaded identify_motif functions, two identify_motif for Hairpins and Internal Loops respectively while identify_motif_b is for bulge loops
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &input_subsequence, char res) {
    if (!initialized){
        initialized = true;
        create_hashmaps();
    }
    Basic_Sequence Motif {&input_subsequence.front(),input_subsequence.size()};
    if (auto search = HairpinHashMap.find(Motif); search != HairpinHashMap.end()){
        return search->second;
    }
    return res;
}

inline char identify_motif(const Basic_Subsequence<char, unsigned int> &internal_subsequence1, const Basic_Subsequence<char, unsigned int> &internal_subsequence2, char res) {
    if (!initialized){
        initialized = true;
        create_hashmaps();
    }   
    Basic_Sequence Motif1 {&internal_subsequence1.front(),internal_subsequence1.size()};
    Basic_Sequence Motif2 {&internal_subsequence2.front(),internal_subsequence2.size()};
    const char connect = char_to_ali_base('$');
    Motif1.append(Motif2.seq,Motif2.size(),&connect);
    if (auto search = InternalHashMap.find(Motif1); search != InternalHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_b(const Basic_Subsequence<char, unsigned int> &bulge_subsequence, char res) {
    if (!initialized){
        initialized = true;
        create_hashmaps();
    }
    Basic_Sequence Motif {&bulge_subsequence.front(),bulge_subsequence.size()};
    if (auto search = BulgeHashMap.find(Motif); search != BulgeHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_align(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, char res) {
    if (!initialized) {
        initialized = true;
        create_hashmaps();
    }
    Basic_Sequence subseq1{&first_track_seq.front(),first_track_seq.size()};
    Basic_Sequence subseq2{&second_track_seq.front(),second_track_seq.size()};
    char found1 = '\0';
    char found2 = '\0';
    char found3;
    char found4;
    if (auto search1 = HairpinHashMap.find(subseq1); search1 != HairpinHashMap.end()) {
        found1 = search1->second;
        if (auto search2 = HairpinHashMap.find(subseq2); search2 != HairpinHashMap.end()) {
            found2 = search2->second;
            if (std::toupper(found1, std::locale()) == std::toupper(found2,std::locale())) {
                return std::toupper(found1, std::locale());
            }
        }
    }
    if (auto search3 = BulgeHashMap.find(subseq1); search3 != BulgeHashMap.end()){
        found3 = search3->second;
        if (auto search4 = BulgeHashMap.find(subseq1); search4 != BulgeHashMap.end()){
            found4 = search4->second;
            if (std::toupper(found3, std::locale()) == std::toupper(found4,std::locale())) {
                return std::toupper(found3, std::locale());
            }
        }
    }
    return res;
}

inline int motif_scoring(const int &length_of_motif_region, const char motif_char) {
    return alignment_match() * length_of_motif_region;
}

//inline bool hairpin_bulge_matching(char FirstFound, char SecondFound){
//    if (std::tolower(FirstFound,std::locale()) == std::tolower(SecondFound,std::locale())){
//        return true;
//    }
//    std::pair values {FirstFound,SecondFound};
//    std::sort(values.first,values.second); //Sort elements, always puts upper case characters first and lower case characters to the back
//    if (std::islower(values.second,std::locale())) { //If the second character isn't lowercase, the first one isn't either.
//        if (std::islower(values.first,std::locale())){
//            std::cout << "bruh\n"; //Double cross check
//        }
//        else {
//            std::cout << "bruh2\n"; //Single cross check
//            return true;
//        }
//    }
//    else{
//        return false;
//    }
//    return true; //dummy to make warning go away
//}
//
//Filter function for RNAmotiFold alignments, makes two hashmap searches and compares the outputs. This ensures only motif matching regions are checked for motifs.
template<typename alphabet, typename pos_type, typename T>
inline bool motif_match(const Basic_Sequence<alphabet, pos_type> &seq1, const Basic_Sequence<alphabet, pos_type> &seq2, T i_seq1, T j_seq1, T i_seq2, T j_seq2){
    Basic_Sequence subseq1{&seq1[i_seq1],j_seq1-i_seq1};
    Basic_Sequence subseq2{&seq2[i_seq2],j_seq2-i_seq2};
    if (!initialized){
        initialized = true;
        create_hashmaps();
    }
    //std::array<bool,3> mot_match;
    //Hairpin Check
    if (auto search = HairpinHashMap.find(subseq1); search != HairpinHashMap.end()){
        char found1 = search->second;
        if (auto search2 = HairpinHashMap.find(subseq2); search2 != HairpinHashMap.end()){
            char found2 = search2->second;
            return std::tolower(found1,std::locale()) == std::tolower(found2,std::locale());
            //mot_match.fill(hairpin_bulge_matching(found1, found2));
            }
    }
    else if (auto search3 = BulgeHashMap.find(subseq1); search3 != BulgeHashMap.end()) {
        char found3 = search3 -> second;
        if (auto search4 = BulgeHashMap.find(subseq2); search4 != BulgeHashMap.end()) {
            char found4 = search4->second;
            return std::tolower(found3,std::locale()) == std::tolower(found4,std::locale());
        }
    }
    else {
        return false;
    }
    return false;
}

//shapeX functions are here to avoid massive if/else statements in the shape_X algebra. This should theoretically make it faster.
//Level [TWO] is always a fallthrough cause there are no motif implementation difference between the two levels.
inline shape_t bl_shapeX(char mot, shape_t &existing_shape){
    const auto level =  static_cast<shapelevel_enum>(gapc::Opts::getOpts()->shapelevel);
    switch (level){
        case five:
            [[fallthrough]];
        case four:
            if (mot != underScore) {
                return shape_t(mot) + existing_shape;
            }
            else{
                return existing_shape;
            }
        case three:
            if (mot != underScore){
                return shape_t(openParen) + shape_t(mot) + existing_shape + shape_t(closeParen);
            }
            else{
                return shape_t(openParen) + existing_shape + shape_t(closeParen);
            }
        case two:
            [[fallthrough]];
        case one:
            return shape_t(openParen) + shape_t(mot) + existing_shape + shape_t(closeParen);
        default:
            std::cerr << "Shape level is not set" << "\n";
            break;
    }
    throw std::invalid_argument("Shape level is not set");
}

inline shape_t br_shapeX(char mot, shape_t &existing_shape){
    const auto level =  static_cast<shapelevel_enum>(gapc::Opts::getOpts()->shapelevel);
    switch (level){
        case five:
            [[fallthrough]];
        case four:
            if (mot != underScore) {
                return existing_shape + shape_t(mot);
            }
            else{
                return existing_shape;
            }
        case three:
            if (mot != underScore){
                return shape_t(openParen) + existing_shape + shape_t(mot) + shape_t(closeParen);
            }
            else{
                return shape_t(openParen) + existing_shape + shape_t(closeParen);
            }
        case two:
            [[fallthrough]];
        case one:
            return shape_t(openParen) + existing_shape + shape_t(mot) + shape_t(closeParen);
        default:
            std::cerr << "Shape level is not set" << "\n";
            break;
    }
    throw std::invalid_argument("Shape level is not set");
}

inline shape_t il_shapeX(char mot, shape_t &existing_shape){
    const auto level =  static_cast<shapelevel_enum>(gapc::Opts::getOpts()->shapelevel);
    switch (level){
        case five:
            if (mot != underScore) {
                return shape_t(mot) + existing_shape + shape_t(mot);
            }
            else{
                return existing_shape;
            }
        case four:
            [[fallthrough]];
        case three:
            if (mot != underScore) {
                return shape_t(openParen) + shape_t(mot) + existing_shape + shape_t(mot) + shape_t(closeParen);
            }
            else {
                return shape_t(openParen) + existing_shape + shape_t(closeParen);
            }
        case two:
            [[fallthrough]];
        case one:
            return shape_t(openParen) + shape_t(mot) + existing_shape + shape_t(mot) +shape_t(closeParen);
        default:
            std::cerr << "Shape level is not set" << "\n";
            break;
    }
    throw std::invalid_argument("Shape level is not set");
}