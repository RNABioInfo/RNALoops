#pragma once
#include "MotifMap.hh"
#include "filter.hh"
#include "motif.hh"
#include "answer_motoh.hh"
#include "rna.hh"
#include "rnaoptions.hh"
#include "rnaoptions_defaults.hh"
#include "sequence.hh"
#include "string.hh"
#include "subsequence.hh"
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include "shape.hh"


namespace motif_ali{
static MotifMap HairpinHashMap {parse_direction_opt()};
static MotifMap InternalHashMap_fronts {parse_direction_opt()};
static MotifMap InternalHashMap_backs {parse_direction_opt()};
static MotifMap InternalHashMap {parse_direction_opt()};
static MotifMap BulgeHashMap {parse_direction_opt()};
static bool init = false;
static std::vector<Basic_Sequence<char, unsigned int>> front_collector_uf1{};
static std::vector<Basic_Sequence<char, unsigned int>> front_collector_uf2{};
static std::vector<Basic_Sequence<char, unsigned int>> front_collector_f1{};
static std::vector<Basic_Sequence<char, unsigned int>> front_collector_f2{};
}

using seq_vector = std::vector<Basic_Subsequence<char, unsigned int>>;

struct MotifsFound {
    std::array <bool,4> states = {false,false,false,false};
    bool& hairpin(){return states[0];};
    bool& bulge(){return states[1];};
    bool& interal_front(){return states[2];};
    bool& internal_full(){return states[3];};
    void setHairpin(bool set_to){states[0] = set_to;};
    void setBulge(bool set_to){states[1] = set_to;};
    void setInternal_front(bool set_to){states[2] = set_to;};
    void setInternal_full(bool set_to){states[3] = set_to;}; 
};

using shape_t = Shape;

template<typename alphabet, typename pos_type>
inline std::set<char> check_maps(Basic_Subsequence<alphabet,pos_type> seq){
    std::set<char> return_set;
    if (auto search = motif_ali::HairpinHashMap.get_motif_set(seq); search != motif_ali::HairpinHashMap.dupe_end()){
        return_set.insert(search->second.begin(),search->second.end());
    }
    if (auto search = motif_ali::InternalHashMap_fronts.get_motif_set(seq); search != motif_ali::InternalHashMap_fronts.dupe_end()){
        return_set.insert(search->second.begin(),search->second.end());
    }
    if (auto search = motif_ali::InternalHashMap_backs.get_motif_set(seq); search != motif_ali::InternalHashMap_backs.dupe_end()){
        return_set.insert(search->second.begin(),search->second.end());
    }
    if (auto search = motif_ali::BulgeHashMap.get_motif_set(seq); search != motif_ali::BulgeHashMap.dupe_end()){
        return_set.insert(search->second.begin(),search->second.end());
    }
    return return_set;
};

template<typename alphabet, typename pos_type>
inline std::set<char> check_map(const Basic_Subsequence<alphabet,pos_type> &seq_first,const Basic_Subsequence<alphabet,pos_type>& seq_second, MotifMap& Checkmap){
    std::set<char> return_set;
    auto search_first = Checkmap.get_motif_set(seq_first);
    auto search_second = Checkmap.get_motif_set(seq_second);
    if (search_first != Checkmap.dupe_end() && search_second != Checkmap.dupe_end()){
        std::set_intersection(search_first->second.begin(),search_first->second.end(),search_second->second.begin(),search_second->second.end(),std::inserter(return_set,return_set.begin()));
    }
    return return_set;
}


inline bool check_for_internal_front(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, const answer_motoh& bruh){
    std::set<char> internal_f = check_map(first_track_seq,second_track_seq,motif_ali::InternalHashMap_fronts);
    if (internal_f.size() > 0){
        return true;
    }
    return false;
}

inline std::set<char> get_motif_overlap(const Basic_Subsequence<char, unsigned int> &seq1, const Basic_Subsequence<char, unsigned int> &seq2){
    std::set<char> Motifs1 = check_maps(seq1);
    std::set<char> Motifs2 = check_maps(seq2);
    std::set<char> res;
    std::set_intersection(Motifs1.begin(),Motifs1.end(),Motifs2.begin(),Motifs2.end(),std::inserter(res,res.begin()));
    return res;
}

inline void append(seq_vector &appended,Basic_Subsequence<char, unsigned int> appendee){
    appended.push_back(appendee);
}

inline char identify_motif_prettier(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq){
    std::set<char> res = get_motif_overlap(first_track_seq,second_track_seq);
    std::set<char>::iterator it = res.begin();
    if (res.size() == 1){
        return *it;
    }
    else{
        return '#';
    }
}

inline String identify_motif_motoh(const Basic_Subsequence<char, unsigned int> &first_track_seq ,const Basic_Subsequence<char, unsigned int> &second_track_seq){
    std::set<char> res = get_motif_overlap(first_track_seq, second_track_seq);
    String return_string;
    return_string.append('<');
    std::set<char>::iterator it = res.begin();
    if (res.size() > 1){
        for (unsigned int index = 0; index < res.size() -1 ; index++){
            return_string.append(*it);
            return_string.append(',');
            std::advance(it,1);
        }
        return_string.append(*it);
    }
    else{
        //return_string.append(*it); Return return string and uncomment this to also include 
        String alternative;
        alternative.empty();
        return alternative;

    }
    return_string.append('>');
    return return_string;
}

inline int identify_motif_mali(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, const answer_motoh& existing_answer) {
    std::set<char> res;
    for (unsigned int index = 0; index < existing_answer.first_track_seqs.size(); index++){
        if (first_track_seq.i - existing_answer.first_track_seqs[index].j  >= 5 && second_track_seq.i - existing_answer.second_track_seqs[index].j >= 5){
        auto search_first = motif_ali::InternalHashMap.get_motif_set(existing_answer.first_track_seqs[index],first_track_seq);
        auto search_second = motif_ali::InternalHashMap.get_motif_set(existing_answer.first_track_seqs[index],second_track_seq);
        if (search_first != motif_ali::InternalHashMap.dupe_end() && search_second != motif_ali::InternalHashMap.dupe_end()) {
            std::set_intersection(search_first->second.begin(),search_first->second.end(),search_second->second.begin(),search_second->second.end(),std::inserter(res,res.begin()));
        }
    }
        else{;}
    }
    if (res.size() > 0) {
        return 10;
        //we found a complete internal loop! let's add a little score
    }
    return 0;
}

inline int motif_scoring(const int &length_of_motif_region) {
    return (alignment_match() * length_of_motif_region) * (alignment_match() * length_of_motif_region);
}

inline std::pair<std::string,std::string> preprocess_internal_motifs(std::string istring){
    constexpr char delim = ',';
    constexpr char internal_delim = '$';
    std::stringstream sinput(istring);
    std::string fronts;
    std::string backs;
    std::string line;
    while (std::getline(sinput,line,'\n')){
        size_t pos = line.find(internal_delim);
        size_t pos2 = line.find(delim);
        fronts.append(line,0,pos);
        fronts.append(line,pos2,std::string::npos);
        fronts.append("\n");
        backs.append(line,pos+1,std::string::npos);
        backs.append("\n");
    }
    std::pair <std::string, std::string> parts{fronts,backs};
    return parts;
}

inline void fill_internal_hashmaps(const std::string & custom_path, bool custom_replace, MotifMap& front_map, MotifMap& back_map, std::array<unsigned char*, 2> arr, std::array<unsigned int,2> len_arr){
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
        motstring.append(add_preset_motifs(arr,len_arr));
    }
    std::pair <std::string,std::string> parts;
    parts = preprocess_internal_motifs(motstring);
    front_map.add_motifs(std::string_view{std::get<0>(parts)});
    back_map.add_motifs(std::string_view{std::get<1>(parts)});
    front_map.print_duplicates();
    back_map.print_duplicates();
}

//Filter function for RNAmotiFold alignments, makes two hashmap searches and compares the outputs. This ensures only motif matching regions are checked for motifs.
template<typename alphabet, typename pos_type, typename T>
inline bool motif_match(const Basic_Sequence<alphabet, pos_type> &seq1, const Basic_Sequence<alphabet, pos_type> &seq2, T i_seq1, T j_seq1, T i_seq2, T j_seq2){
    if (!motif_ali::init) {
        fill_hashmap(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts()->replaceH, motif_ali::HairpinHashMap, motif_basic::Hairpins, motif_basic::Hairpin_lengths);
        fill_hashmap(gapc::Opts::getOpts()->custom_bulges, gapc::Opts::getOpts()->replaceB, motif_ali::BulgeHashMap, motif_basic::Bulges, motif_basic::Bulge_lengths);
        fill_hashmap(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts()->replaceI, motif_ali::InternalHashMap, motif_basic::Internals,motif_basic::Internal_lengths);
        fill_internal_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts()->replaceI, motif_ali::InternalHashMap_fronts, motif_ali::InternalHashMap_backs, motif_basic::Internals, motif_basic::Internal_lengths);
        motif_ali::init = true;
    }
    Basic_Subsequence<alphabet, pos_type> Seq1 {seq1, i_seq1, j_seq1};
    Basic_Subsequence<alphabet, pos_type> Seq2 {seq2, i_seq2, j_seq2};
    std::set<char> res;
    if (j_seq1 < seq1.size() -1 && j_seq2 < seq2.size() -1 && i_seq1 > 1 && i_seq2 > 1){ //Bei j size -1 damit theoretisch zwei basepairs noch hinpassen und bei i größer 1 damit theoretisch noch 2 basepairs davor passen
        if (rnali_basepairing(seq1, i_seq1 - 1, j_seq1 + 1) && rnali_basepairing(seq2, i_seq2 - 1, j_seq2 + 1) && rnali_basepairing(seq1, i_seq1 - 2 , j_seq1 + 2) && rnali_basepairing(seq2, i_seq2 - 2, j_seq2 + 2)){ // Sowohl i+1/j+1 und i+2/j+2 müssen pairen? Mit Lonely Basepairs ?
            res.merge(check_map(Seq1, Seq2 , motif_ali::HairpinHashMap));
        }
    }
    std::set<char> internal_fronts = check_map(Seq1,Seq2, motif_ali::InternalHashMap_fronts);
    if (internal_fronts.size() > 0){
        res.merge(internal_fronts);
        Basic_Sequence Motif1 {&Seq1.front(),Seq1.size()};
        Basic_Sequence Motif2 {&Seq2.front(),Seq2.size()};
        motif_ali::front_collector_uf1.push_back(Motif1);
        motif_ali::front_collector_uf2.push_back(Motif2);
    }
    std::set<char> internal_backs = check_map(Seq1,Seq2, motif_ali::InternalHashMap_backs);
    if (internal_backs.size() > 0){
        res.merge(internal_backs);
        Basic_Sequence Motif3 {&Seq1.front(),Seq1.size()};
        Basic_Sequence Motif4 {&Seq2.front(),Seq2.size()};
        for (auto seq: motif_ali::front_collector_uf1){
            auto cpy = seq;
            cpy.concat(Motif3, Motif3.size());
            if (auto search = motif_ali::InternalHashMap.get_motif_set(cpy); search != motif_ali::InternalHashMap.dupe_end()){
                motif_ali::front_collector_f1.push_back(seq);
            }
        }
    }

    return res.size() > 0;
    }
