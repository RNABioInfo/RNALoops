#pragma once
#include "subsequence.hh"
#include "sequence.hh"
#include "rnaoptions.hh"
#include "mot_header.hh"
#include "shapes.hh"
#include "shape.hh"
#include <cctype>
#include <iterator>
#include <stdexcept>
#include <string>
#include <iostream>
#include <fstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include"MotifMap.hh"

struct HairpinLoopMotif{};
struct InternalLoopMotif{};
struct BulgeLoopMotif{};

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
    bool initializedH(){return states[0];};
    bool initializedI(){return states[1];};
    bool initializedB(){return states[2];};
    void setH(bool set){states[0] = set;};
    void setI(bool set){states[1] = set;};  
    void setB(bool set){states[2] = set;};
    void setAll(bool set){setH(set);setI(set);setB(set);};
    bool initialized() {return std::all_of(states.begin(),states.end(),[](const bool v){return v;});};
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

template <typename T>
inline void initializer(T loop_type){
    if (!init.initialized()){
        initialize_hash_map(init, loop_type);
    }
}

template<typename T>
inline void initialize_hash_map(init_status &init_check, T _){
    if constexpr(std::is_same<HairpinLoopMotif, T>::value){
        if (!init_check.initializedH()){
            create_hashmaps(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts() -> replaceH, HairpinHashMap, Hairpins, Hairpin_lengths);
            init_check.setH(true);
        }
        return;
    }
    else if constexpr(std::is_same<InternalLoopMotif, T>::value){
        if (!init_check.initializedI()){
            create_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, InternalHashMap, Internals, Internal_lengths);
            init_check.setI(true);
        }
        return; 
    }
    else if constexpr(std::is_same<BulgeLoopMotif, T>::value){
        if (!init_check.initializedB()){
            create_hashmaps(gapc::Opts::getOpts()->custom_bulges, gapc::Opts::getOpts() -> replaceB, BulgeHashMap, Bulges, Bulge_lengths);
            init_check.setB(true);
        }
        return;
    }
    else {
        throw std::runtime_error("Unknown Motif Type during Hash Map Initialization");
    }
}

template<typename pos_type>
inline char select_return_motif(std::vector<char> found_motifs, pos_type rows , char res){
    CounterMap Map(found_motifs);
    std::pair<char, unsigned int> maxval = Map.findMaxValuePair();
    double fract = gapc::Opts::getOpts()->fraction;
    if (static_cast<double>(maxval.second)/rows >= fract) {
         return maxval.first;
    }
    else{
        return res;
    }
}

template<typename pos_type>
inline int select_return_motif(std::vector<char> found_motifs, pos_type rows, std::unordered_map<char, std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array>> seq_versions){
    CounterMap Map(found_motifs);
    std::pair<char, unsigned int> maxval = Map.findMaxValuePair();
    double fract = gapc::Opts::getOpts()->fraction;
    float weight = gapc::Opts::getOpts()->weighting;
    if (static_cast<double>(maxval.second)/rows >= fract) {
            double versions = static_cast<double>(seq_versions[maxval.first].size());
            if (versions == 1){
                return 1.0 * weight;
            }
            else{
                return -std::ceil((versions/static_cast<double>(rows))*weight);
            }
        }
}

inline std::vector<char> find_and_count_motifs(std::unordered_map<char, std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array>> &versions, const Basic_Subsequence<M_Char, unsigned int> &seq, MotifMap &Map){
    std::vector<char> found;
    for (unsigned row = 0; row < rows(seq); row++){
        std::vector<long unsigned int> vec = get_gaps(seq.seq->row(row),seq.i,seq.j);
        Basic_Sequence Motif{seq.seq->row(row), seq.i,seq.j,vec};
        if (auto search = Map.find(Motif); search != Map.end()) {
            for (unsigned int i = 0; i < Map.Dupes[Motif].size(); i++){
                char mots = *next(Map.Dupes[Motif].begin(),i);
                found.push_back(mots);
                if (versions.find(mots) == versions.end()){
                    versions[mots] = std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array> {Motif};
                }
                else{
                    versions[mots].insert(Motif);
                }
            }
        }
    }
    return found;
}

inline std::vector<char> find_and_count_motifs(std::unordered_map<char, std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array>> &versions, const Basic_Subsequence<M_Char, unsigned int> &seq, const Basic_Subsequence<M_Char, unsigned int> &seq2, MotifMap &Map){
    std::vector<char> found;
    for (unsigned row = 0; row < rows(seq); row++){
            std::vector<long unsigned int> vec = get_gaps(seq.seq->row(row),seq.i,seq.j);
            std::vector<long unsigned int> vec2 = get_gaps(seq2.seq->row(row),seq2.i,seq2.j);
            Basic_Sequence Motif{seq.seq->row(row), seq.i,seq.j,vec};
            Basic_Sequence Motif2{seq2.seq->row(row), seq2.i,seq2.j,vec2};
            Motif.concat(Motif2.seq,Motif2.size());
        if (auto search = Map.find(Motif); search != Map.end()) {
            for (unsigned int i = 0; i < Map.Dupes[Motif].size(); i++){
                char mots = *next(Map.Dupes[Motif].begin(),i);
                found.push_back(mots);
                if (versions.find(mots) == versions.end()){
                    versions[mots] = std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array> {Motif};
                }
                else{
                    versions[mots].insert(Motif);
                }
            }
        }
    }
    return found;
}

//Regular identify motif functions for either Single Sequence Folding with Basic_Subsequences or Alignment Folding with M_Char
template<typename T>
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &motif_sequence, char res, T loop_type){
    initializer(loop_type);
    if constexpr(std::is_same_v<HairpinLoopMotif, T>){
        if (auto search = HairpinHashMap.find(motif_sequence); search != HairpinHashMap.end()){
            return search->second;
        }
    }
    else if constexpr(std::is_same_v<BulgeLoopMotif, T>){
        if (auto search = BulgeHashMap.find(motif_sequence); search != BulgeHashMap.end()) {
            return search->second;  
        }
    }
    else{
        throw std::runtime_error("Unknown Motif Type detected during motif identification!");
    }
    return res;
}

template<typename T>
inline char identify_motif(const Basic_Subsequence<char, unsigned int> &internal_subsequence1, const Basic_Subsequence<char, unsigned int> &internal_subsequence2, char res, T loop_type){
    initializer(loop_type);
    if (auto search = InternalHashMap.find(internal_subsequence1,internal_subsequence2); search != InternalHashMap.end()) {
        return search->second;  
    }
    return res;
}

template <typename T>
inline char identify_motif(const Basic_Subsequence<M_Char, unsigned int> &motif_sequence,char res, T loop_type) {
    initializer(loop_type);
    std::vector<char> found;
    if constexpr(std::is_same_v<HairpinLoopMotif, T>){
        found = HairpinHashMap.find(motif_sequence);
    }
    else if constexpr(std::is_same_v<BulgeLoopMotif, T>){
        found = BulgeHashMap.find(motif_sequence);
    }
  if (!found.empty()) {
    return select_return_motif(found, rows(motif_sequence), res);
  }
  return res;
}

template <typename T>
inline char identify_motif(const Basic_Subsequence<M_Char, unsigned int> &internal_subsequence1,const Basic_Subsequence<M_Char, unsigned int> &internal_subsequence2,char res, T loop_type) {
    initializer(loop_type);
    std::vector<char> found;
    found = InternalHashMap.find(internal_subsequence1,internal_subsequence2);
    if (!found.empty()) {
        return select_return_motif(found, rows(internal_subsequence1), res);
    }
    return res;
}

//Scoring function for RNAmotiAlign, returns Score for a set of Basic_Subsequence objects
template <typename T>
inline float motifscore(const Basic_Subsequence<M_Char, unsigned int> &seq, T loop_type){
    initializer(loop_type);
    std::vector<char> found;
    std::unordered_map<char, std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array>> seq_versions{};
    if constexpr (std::is_same_v<HairpinLoopMotif,T>){
        found = find_and_count_motifs(seq_versions, seq,HairpinHashMap);
    }
    else if constexpr(std::is_same_v<BulgeLoopMotif,T>){
        found = find_and_count_motifs(seq_versions, seq,BulgeHashMap);
    }
    if (!found.empty()){        
        return select_return_motif(found, rows(seq),seq_versions);
    }
    return 0;
}

template <typename T>
inline float motifscore(const Basic_Subsequence<M_Char, unsigned int> &seq,const Basic_Subsequence<M_Char, unsigned int> &seq2, T loop_type){
    initializer(loop_type);
    std::vector<char> found;
    std::unordered_map<char, std::unordered_set<Basic_Sequence<char,unsigned int>,Hash_ali_array>> seq_versions{};
    found = find_and_count_motifs(seq_versions, seq,seq2,InternalHashMap);
    if (!found.empty()){        
        return select_return_motif(found, rows(seq),seq_versions);
    }
    return 0;
}

//Syntactic filter function for hairpin loop motifs (with(_overlay), pre-parsing)
template <typename alphabet, typename pos_type, typename T>
inline bool motif_h(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
    if (!init.initializedH()){
        create_hashmaps(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts() -> replaceH, HairpinHashMap, Hairpins, Hairpin_lengths);
        init.setH(true);
    }
    if (auto search = HairpinHashMap.find(seq,i,j); search != HairpinHashMap.end()) {
    return true;
  }
  return false;
}

//Syntactic filter function for bulge loop motifs (with(_overlay), pre-parsing)
template <typename alphabet, typename pos_type, typename T>
inline bool motif_b(const Basic_Sequence<alphabet, pos_type> &seq, T i, T j) {
    if (!init.initializedB()){
        create_hashmaps(gapc::Opts::getOpts()->custom_bulges, gapc::Opts::getOpts() -> replaceB, BulgeHashMap, Bulges, Bulge_lengths);
        init.setB(true);
    }
    if (auto search = BulgeHashMap.find(seq,i,j); search != BulgeHashMap.end()) {
        return true;
    }
    return false;
}

//Syntactic overlay filter function for internal loop motifs (with(_overlay) , pre-parsing)
template<typename alphabet, typename pos_type, typename T>
inline bool motif_i (const  Basic_Sequence<alphabet, pos_type> &seq, T lb_i, T lb_j, T lr_i, T lr_j, T x_i, T x_j, T rr_i, T rr_j, T rb_i, T rb_j){
    if (!init.initializedI()){
        create_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, InternalHashMap, Internals, Internal_lengths);
        init.setI(true);
    }
    const Basic_Subsequence<alphabet,pos_type>& lr {seq,lr_i,lr_j};
    const Basic_Subsequence<alphabet,pos_type>& rr {seq,rr_i,rr_j};
    if (auto search = InternalHashMap.find(lr,rr); search != InternalHashMap.end()){
        return true;
    }
    return false;
}

//Semantic overlay filter function for motif_i (suchthat(_overlay), post-parsing). Functional but currently not in use cause this filter is algebra independent!
template <typename alphabet, typename pos_type, typename T>
inline bool motif_i(const Basic_Subsequence<alphabet, pos_type> &base1,
                    const Basic_Subsequence<alphabet, pos_type> &seq1, 
                    T & _ ,
                    const Basic_Subsequence<alphabet, pos_type> &seq2,
                    const Basic_Subsequence<alphabet, pos_type> &base2) {
    if (!init.initializedI()){
        create_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts() -> replaceI, InternalHashMap, Internals, Internal_lengths);
        init.setI(true);
    }
    if (auto search = InternalHashMap.find(seq1, seq2); search != InternalHashMap.end()) {
        return true;
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