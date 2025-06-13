#pragma once
#include "motif.hh"
#include "subsequence.hh"
#include <algorithm>
#include <stdexcept>
#include "shape.hh"


namespace motif_ali{
static MotifMap HairpinHashMap {parse_direction_opt()};
static MotifMap InternalHashMap1 {parse_direction_opt()};
static MotifMap InternalHashMap2 {parse_direction_opt()};
static MotifMap BulgeHashMap {parse_direction_opt()};
static bool ali_init = false;
}

using shape_t = Shape;

inline char identify_motif_motoh(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq) {
    if (auto search1 = motif_ali::AlignmentMap.get_motif(first_track_seq); search1 != AlignmentMap.end()) {
        char found1 = search1->second;
        if (auto search2 = AlignmentMap.get_motif(second_track_seq); search2 != AlignmentMap.end()) {
            char found2 = search2->second;
            if (std::tolower(found1,std::locale())== std::tolower(found2,std::locale())){
                return found1;
                }
            }
        }
    throw std::runtime_error("Previously identified motif region could not be identified again ? How did we get here ?");
}

inline std::pair<std::pair<Basic_Subsequence<char, unsigned int>,Basic_Subsequence<char, unsigned int>>,std::set<char>> identify_motif_mali(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq) {
    std::set<char> res;
    std::vector<std::set<char>> collect;
    if (auto search1 = AlignmentMap.get_motif_set(first_track_seq); search1 != AlignmentMap.dupe_end()) {
        std::set <char> found1 = search1->second;
        if (auto search2 = AlignmentMap.get_motif_set(second_track_seq); search2 != AlignmentMap.dupe_end()) {
            std::set <char> found2 = search2->second;  
            std::set_intersection(found1.begin(),found1.end(),found2.begin(),found2.end(),std::inserter(res,res.begin()));
        }
    }
    return std::pair<std::pair<Basic_Subsequence<char, unsigned int>,Basic_Subsequence<char, unsigned int>>,std::set<char>> {std::pair<Basic_Subsequence<char, unsigned int>,Basic_Subsequence<char, unsigned int>>{first_track_seq,second_track_seq},res};
}

inline int motif_scoring(const int &length_of_motif_region) {
    return (alignment_match() * length_of_motif_region) + 1;
}

//Filter function for RNAmotiFold alignments, makes two hashmap searches and compares the outputs. This ensures only motif matching regions are checked for motifs.
template<typename alphabet, typename pos_type, typename T>
inline bool motif_match(const Basic_Sequence<alphabet, pos_type> &seq1, const Basic_Sequence<alphabet, pos_type> &seq2, T i_seq1, T j_seq1, T i_seq2, T j_seq2){
    if (!motif_ali::ali_init) {
        fill_hashmap(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts()->replaceH, motif_ali::HairpinHashMap, motif_basic::Hairpins, motif_basic::Hairpin_lengths);
        fill_hashmap(gapc::Opts::getOpts()->custom_bulges, gapc::Opts::getOpts()->replaceB, motif_ali::BulgeHashMap, motif_basic::Bulges, motif_basic::Bulge_lengths);
        motif_ali::ali_init = true;
    }
    Basic_Subsequence<char, unsigned int> Motif1 {seq1,i_seq1,j_seq1};
    Basic_Subsequence<char, unsigned int> Motif2 {seq2,i_seq2,j_seq2};
    if (auto search = AlignmentMap.get_motif_set(Motif1); search != AlignmentMap.dupe_end()){
        std::set<char> found1 = search->second;
        if (auto search2 = AlignmentMap.get_motif_set(Motif2); search2 != AlignmentMap.dupe_end()){
            std::set<char> found2 = search2->second;
            //set comparison function that checks if at least one of the motifs overlaps:
            std::set<char> res;
            std::set_intersection(found1.begin(),found1.end(),found2.begin(),found2.end(),std::inserter(res,res.begin()));
            return res.size() > 0;
            }
    }
    else {
        return false;
    }
    return false;
}