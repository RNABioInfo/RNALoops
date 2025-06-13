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

inline static  direction_type parse_direction_opt(){
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


using shape_t =  Shape;
namespace motif_basic{
static MotifMap HairpinHashMap{parse_direction_opt()};
static MotifMap InternalHashMap{parse_direction_opt()};
static MotifMap BulgeHashMap{parse_direction_opt()};
static std::array Hairpins         = {rna3d_hairpins, rfam_hairpins};
static std::array Hairpin_lengths  = {rna3d_hairpins_len, rfam_hairpins_len};
static std::array Internals        = {rna3d_internals, rfam_internals};
static std::array Internal_lengths = {rna3d_internals_len, rfam_internals_len,};
static std::array Bulges           = {rna3d_bulges, rfam_bulges};
static std::array Bulge_lengths    = {rna3d_bulges_len, rfam_bulges_len};
}
enum shapelevel_enum: std::uint8_t {five=5,four=4,three=3,two=2,one=1};

struct init_status {
    std::array <bool,3> states = {false,false,false};
    bool& initializedH(){return states[0];};
    bool& initializedI(){return states[1];};
    bool& initializedB(){return states[2];};
    void setHairpin(bool set_to){states[0] = set_to;};
    void setInternal(bool set_to){states[1] = set_to;};  
    void setBulge(bool set_to){states[2] = set_to;};
    void setAll(bool set_all_to){setHairpin(set_all_to);setInternal(set_all_to);setBulge(set_all_to);};
};
static init_status init;

//HashMap implementation functions for all three loop types in the macrostate grammar.DB =database, which one gets used 1=BGSU, 2=RMFAM, 3 = both // R = reverse, 1= No Reverses, 2 = Only Reverses, 3 = Both reverse and Forward
inline std::string get_preset_motifs(std::array<unsigned char* ,2> motif_array, std::array<unsigned int,2> motif_len_array) {
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
    return str_mot_ar;
}

//function that alters the global HashMaps (I know this might be bad practice but it's easy and it works.)
inline void fill_hashmap(const std::string& custom_path, bool custom_replace, MotifMap& empty_map ,std::array<unsigned char* ,2> arr, std::array<unsigned int,2> len_arr){
    std::string motstring;
    if (!custom_path.empty()){
        std::ifstream ifstream(custom_path);
        if (!ifstream.is_open()){
            throw std::runtime_error("File not found");
        }
        motstring.append((std::istreambuf_iterator<char>(ifstream)),std::istreambuf_iterator<char>());
        motstring.append("\n");
    }
    if (!custom_replace){
        motstring.append(get_preset_motifs(arr,len_arr));
    }
    std::string_view mot_view {motstring};
    empty_map.add_motifs(mot_view);
    empty_map.print_duplicates();
}

//Overloaded identify_motif functions, two identify_motif for Hairpins and Internal Loops respectively while identify_motif_b is for bulge loops
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &input_subsequence, char res) {
    if (!init.initializedH()){
        init.setHairpin(true);
        fill_hashmap(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts() -> replaceH, motif_basic::HairpinHashMap,motif_basic::Hairpins,motif_basic::Hairpin_lengths);
    }
    if (auto search = motif_basic::HairpinHashMap.get_motif(input_subsequence); search != motif_basic::HairpinHashMap.end()){
        return search->second;
    }
    return res;
}

inline char identify_motif(const Basic_Subsequence<char, unsigned int> &internal_subsequence1, const Basic_Subsequence<char, unsigned int> &internal_subsequence2, char res) {
    if (!init.initializedI()){
        init.setInternal(true);
        fill_hashmap(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, motif_basic::InternalHashMap,motif_basic::Internals,motif_basic::Internal_lengths);
    }
    if (auto search = motif_basic::InternalHashMap.get_motif(internal_subsequence1,internal_subsequence2); search != motif_basic::InternalHashMap.end()) {
        return search->second;
    }
    return res;
}

inline char identify_motif_b(const Basic_Subsequence<char, unsigned int> &bulge_subsequence, char res) {
    if (!init.initializedB()){
        init.setBulge(true);
        fill_hashmap(gapc::Opts::getOpts() -> custom_bulges,gapc::Opts::getOpts() -> replaceB, motif_basic::BulgeHashMap,motif_basic::Bulges,motif_basic::Bulge_lengths);
    }
    if (auto search = motif_basic::BulgeHashMap.get_motif(bulge_subsequence); search != motif_basic::BulgeHashMap.end()) {
        return search->second;
    }
    return res;
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