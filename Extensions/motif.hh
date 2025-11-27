#pragma once
#include "rnaoptions_defaults.hh"
#include "subsequence.hh"
#include "sequence.hh"
#include "rnaoptions.hh"
#include "mot_header.hh"
#include "shapes.hh"
#include "shape.hh"
#include <cctype>
#include <iterator>
#include <locale>
#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <string_view>
#include <utility>
#include <algorithm>
#include"MotifMap.hh"


using shape_t =  Shape;
static MotifMap HairpinHashMap;
static MotifMap InternalHashMap;
static MotifMap BulgeHashMap;
static std::array Hairpins         = {rna3d_hairpins, rfam_hairpins};
static std::array Hairpin_lengths  = {rna3d_hairpins_len, rfam_hairpins_len};
static std::array Internals        = {rna3d_internals, rfam_internals};
static std::array Internal_lengths = {rna3d_internals_len, rfam_internals_len,};
static std::array Bulges           = {rna3d_bulges, rfam_bulges};
static std::array Bulge_lengths    = {rna3d_bulges_len, rfam_bulges_len};

enum shapelevel_enum: std::uint8_t {five=5,four=4,three=3,two=2,one=1};
struct init_status {
    std::array <bool,3> states = {false,false,false};
    bool& initializedH(){return states[0];};
    bool& initializedI(){return states[1];};
    bool& initializedB(){return states[2];};
    void setH(bool set){states[0] = set;};
    void setI(bool set){states[1] = set;};  
    void setB(bool set){states[2] = set;};
    void setAll(bool set){setH(set);setI(set);setB(set);};
};
static init_status init;


inline static direction_type parse_direction_opt(){
    switch(gapc::Opts::getOpts()->reversed){
        case 1:
            return direction_type::forward;
        case 2:
            return direction_type::reverse;
        case 3:
            return direction_type::both;
        default:
            throw std::runtime_error("Something went wrong with the directions...");
    }
};

//HashMap implementation functions for all three loop types in the macrostate grammar.DB =database, which one gets used 1=BGSU, 2=RMFAM, 3 = both // R = reverse, 1= No Reverses, 2 = Only Reverses, 3 = Both reverse and Forward
inline MotifMap Motif_HashMap(std::array<unsigned char* ,2> motif_array, std::array<unsigned int,2> motif_len_array, direction_type direction) {
    std::string str_mot_ar;
    switch(gapc::Opts::getOpts()->motifs) {
        case 1:
            str_mot_ar.append(reinterpret_cast<char*>(motif_array[0]),motif_len_array[0]);
            break;
        case 2:
            str_mot_ar.append(reinterpret_cast<char*>(motif_array[1]),motif_len_array[1]);
            break;
        case 3:
            str_mot_ar.append(reinterpret_cast<char*>(motif_array[0]),motif_len_array[0]);
            str_mot_ar.append(reinterpret_cast<char*>(motif_array[1]),motif_len_array[1]);
            break;
        default:
            throw std::runtime_error("Could not identify given built in motif sets");
    }
    std::string_view istring {str_mot_ar};
    MotifMap NewMap (istring,direction);
    NewMap.print_duplicates();
    return NewMap;
};

//function that alters the global HashMaps (I know this might be bad practice but it's easy and it works.)
inline void fill_hashmap(const std::string& custom_path, bool custom_replace, MotifMap& empty_map ,std::array<unsigned char* ,2> arr, std::array<unsigned int,2> len_arr, direction_type directions){
    if (!custom_path.empty()){
        std::ifstream ifstream(custom_path);
        if (!ifstream.is_open()){
            throw std::runtime_error("File not found");
        }
        std::string custom_motifs((std::istreambuf_iterator<char>(ifstream)),std::istreambuf_iterator<char>());
        std::string_view istring {custom_motifs};

        if (custom_replace){
            empty_map.add_motifs(istring);
        }
        else{
            empty_map = Motif_HashMap(arr,len_arr,directions);
            empty_map.add_motifs(istring);
        }
    }
    else{
        empty_map = Motif_HashMap(arr,len_arr,directions);
    }
}

inline void create_hashmaps(const std::string& custom_path, bool replace_bool, MotifMap& map_to_fill, std::array<unsigned char*,2> motif_arrays, std::array<unsigned int,2> length_arrays){
    direction_type direct = parse_direction_opt();
    fill_hashmap(custom_path,    replace_bool,    map_to_fill, motif_arrays,    length_arrays, direct); 
}

inline void initializeH(init_status& init_check){
    if (!init_check.initializedH()){
        init_check.setH(true);
        create_hashmaps(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts() -> replaceH, HairpinHashMap, Hairpins, Hairpin_lengths);
    }
};

inline void initializeI(init_status& init_check){
    if (!init_check.initializedI()){
        init_check.setI(true);
        create_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, InternalHashMap, Internals, Internal_lengths);
    }
};

inline void initializeB(init_status& init_check){
    if (!init_check.initializedB()){
        init_check.setB(true);
        create_hashmaps(gapc::Opts::getOpts()->custom_bulges, gapc::Opts::getOpts() -> replaceB, BulgeHashMap, Bulges, Bulge_lengths);
    }
};

//Overloaded identify_motif functions, two identify_motif for Hairpins and Internal Loops respectively while identify_motif_b is for bulge loops
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &input_subsequence, char res) {
    initializeH(init);
    if (auto search = HairpinHashMap.find(input_subsequence); search != HairpinHashMap.end()){
        return search->second;
    }
    return res;
}

inline char identify_motif(const Basic_Subsequence<char, unsigned int> &internal_subsequence1, const Basic_Subsequence<char, unsigned int> &internal_subsequence2, char res) {
    initializeI(init);
    if (auto search = InternalHashMap.find(internal_subsequence1,internal_subsequence2); search != InternalHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_b(const Basic_Subsequence<char, unsigned int> &bulge_subsequence, char res) {
    initializeB(init);
    if (auto search = BulgeHashMap.find(bulge_subsequence); search != BulgeHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_align(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, char res) {
    if (!std::all_of(init.states.cbegin(),init.states.cend(),[](bool init){return init;})) {
        create_hashmaps(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts() -> replaceH, HairpinHashMap,Hairpins,Hairpin_lengths);
        create_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, InternalHashMap,Internals,Internal_lengths);
        create_hashmaps(gapc::Opts::getOpts() -> custom_bulges,gapc::Opts::getOpts() -> replaceB, BulgeHashMap,Bulges,Bulge_lengths);
        init.setAll(true);
    }
    char found1;
    char found2;
    char found3;
    char found4;
    if (auto search1 = HairpinHashMap.find(first_track_seq); search1 != HairpinHashMap.end()) {
        found1 = search1->second;
        if (auto search2 = HairpinHashMap.find(second_track_seq); search2 != HairpinHashMap.end()) {
            found2 = search2->second;
            if (std::toupper(found1, std::locale()) == std::toupper(found2,std::locale())) {
                return std::toupper(found1, std::locale());
            }
        }
    }
    if (auto search3 = BulgeHashMap.find(first_track_seq); search3 != BulgeHashMap.end()){
        found3 = search3->second;
        if (auto search4 = BulgeHashMap.find(second_track_seq); search4 != BulgeHashMap.end()){
            found4 = search4->second;
            if (std::toupper(found3, std::locale()) == std::toupper(found4,std::locale())) {
                return std::toupper(found3, std::locale());
            }
        }
    }
    return res;
}

//Syntactic filter function for hairpin loop motifs (with(_overlay), pre-parsing)
template <typename alphabet, typename pos_type, typename T>
inline bool motif_h(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
    initializeH(init);
    if (auto search = HairpinHashMap.find(seq,i,j); search != HairpinHashMap.end()) {
    return true;
  }
  return false;
}

//Syntactic filter function for bulge loop motifs (with(_overlay), pre-parsing)
template <typename alphabet, typename pos_type, typename T>
inline bool motif_b(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
    initializeB(init);
    if (auto search = BulgeHashMap.find(seq,i,j); search != BulgeHashMap.end()) {
        return true;
    }
    return false;
}

//Semantic overlay filter function for motif_i (suchthat(_overlay), post-parsing). Functional but currently not in use cause this filter is algebra independent!
template <typename alphabet, typename pos_type, typename T>
inline bool motif_i(const Basic_Subsequence<alphabet, pos_type> &base1,
                    const Basic_Subsequence<alphabet, pos_type> &seq1, 
                    T &bruh,
                    const Basic_Subsequence<alphabet, pos_type> &seq2,
                    const Basic_Subsequence<alphabet, pos_type> &base2) {
    initializeI(init);
    if (auto isearch = InternalHashMap.find(seq1, seq2); isearch != InternalHashMap.end()) {
        return true;
    }
    return false;
}

//Syntactic overlay filter function for internal loop motifs (with(_overlay) , pre-parsing)
template<typename alphabet, typename pos_type, typename T>
inline bool motif_i (const  Basic_Sequence<alphabet, pos_type> &seq, T lb_i, T lb_j, T lr_i, T lr_j, T x_i, T x_j, T rr_i, T rr_j, T rb_i, T rb_j){
    initializeI(init);
    const Basic_Subsequence<alphabet,pos_type>& lr {seq,lr_i,lr_j};
    const Basic_Subsequence<alphabet,pos_type>& rr {seq,rr_i,rr_j};
    if (auto search = InternalHashMap.find(lr,rr); search != InternalHashMap.end()){
        return true;
    }
    return false;
}

//Filter function for RNAmotiFold alignments, makes two hashmap searches and compares the outputs. This ensures only motif matching regions are checked for motifs.
template<typename alphabet, typename pos_type, typename T>
inline bool motif_match(const Basic_Sequence<alphabet, pos_type> &seq1, const Basic_Sequence<alphabet, pos_type> &seq2, T i_seq1, T j_seq1, T i_seq2, T j_seq2){
    if (!std::all_of(init.states.cbegin(),init.states.cend(),[](bool init){return init;})) {
        create_hashmaps(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts() -> replaceH, HairpinHashMap,Hairpins,Hairpin_lengths);
        create_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, InternalHashMap,Internals,Internal_lengths);
        create_hashmaps(gapc::Opts::getOpts() -> custom_bulges,gapc::Opts::getOpts() -> replaceB, BulgeHashMap,Bulges,Bulge_lengths);
        init.setAll(true);
    }
    //Hairpin Check
    Basic_Subsequence<char, unsigned int> Motif1 {seq1,i_seq1,j_seq1};
    Basic_Subsequence<char, unsigned int> Motif2 {seq2,i_seq2,j_seq2};
    if (auto search = HairpinHashMap.find(Motif1); search != HairpinHashMap.end()){
        char found1 = search->second;
        if (auto search2 = HairpinHashMap.find(Motif2); search2 != HairpinHashMap.end()){
            char found2 = search2->second;
            return std::tolower(found1,std::locale()) == std::tolower(found2,std::locale());
            }
    }
    else if (auto search3 = BulgeHashMap.find(Motif1); search3 != BulgeHashMap.end()) {
        char found3 = search3 -> second;
        if (auto search4 = BulgeHashMap.find(Motif2); search4 != BulgeHashMap.end()) {
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