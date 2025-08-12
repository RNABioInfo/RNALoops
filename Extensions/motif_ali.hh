#pragma once
#include "MotifMap.hh"
#include "filter.hh"
#include "motif.hh"
#include "backtrack.hh"
#include "answer_motoh.hh"
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
#include "shape.hh"

using Internal_location = std::pair<std::pair<unsigned int, unsigned int>, std::pair <unsigned int, unsigned int>>;
namespace motif_ali{
static MotifMap HairpinHashMap {parse_direction_opt()};
static MotifMap InternalHashMap_fronts {parse_direction_opt()};
static MotifMap InternalHashMap_backs {parse_direction_opt()};
static MotifMap InternalHashMap {parse_direction_opt()};
static MotifMap BulgeHashMap {parse_direction_opt()};
static bool init = false;
static std::vector<Internal_location> back_collector{};
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

inline bool identify_internal_back(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, const answer_motoh& bruh){
    std::set<char> internal_b = check_map(first_track_seq,second_track_seq,motif_ali::InternalHashMap_backs);
    return internal_b.size() > 0;
}
inline bool identify_internal_front(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, const answer_motoh& bruh){
    std::set<char> internal_f = check_map(first_track_seq,second_track_seq,motif_ali::InternalHashMap_fronts);
    return internal_f.size() > 0;
}

inline std::set<char> get_motif_overlap(const Basic_Subsequence<char, unsigned int> &seq1, const Basic_Subsequence<char, unsigned int> &seq2){
    std::set<char> Motifs1 = check_maps(seq1);
    std::set<char> Motifs2 = check_maps(seq2);
    std::set<char> res;
    std::set_intersection(Motifs1.begin(),Motifs1.end(),Motifs2.begin(),Motifs2.end(),std::inserter(res,res.begin()));
    return res;
}

inline void append(seq_vector &appended,const Basic_Subsequence<char, unsigned int> appendee){
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
    std::set<char>::iterator it = res.begin();
    if (res.size() > 1){
        return_string.append('<');
        for (unsigned int index = 0; index < res.size() -1 ; index++){
            return_string.append(*it);
            return_string.append(',');
            std::advance(it,1);
        }
        return_string.append(*it);
        return_string.append('>');
    }
    else{
        return_string.append(*it);// Return return string and uncomment this to also include single letters ig ?
        //String alternative;
        //alternative.empty();
        //return alternative;
    }
    
    return return_string;
}

inline int identify_motif_mali(const Basic_Subsequence<char, unsigned int> &first_track_seq, const Basic_Subsequence<char, unsigned int> &second_track_seq, const answer_motoh& existing_answer) {
    std::set<char> res;
    for (unsigned int index = 0; index < existing_answer.first_track_seqs.size(); index++){
            auto search_first = motif_ali::InternalHashMap.get_motif_set(first_track_seq,existing_answer.first_track_seqs[index]); //Combine the new found front with the previously found back for Seq1
            auto search_second = motif_ali::InternalHashMap.get_motif_set(second_track_seq,existing_answer.second_track_seqs[index]); //Combine the new found front with the previously found back for Seq2
            if (search_first != motif_ali::InternalHashMap.dupe_end() && search_second != motif_ali::InternalHashMap.dupe_end()) {
                std::set_intersection(search_first->second.begin(),search_first->second.end(),search_second->second.begin(),search_second->second.end(),std::inserter(res,res.begin()));
        }
    }
    return res.size() > 0;
}

inline int motif_scoring(const int &length_of_motif_region) {
    return (alignment_match() * length_of_motif_region);
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

template<typename T>
struct equilibrium {
    int val;
    equilibrium(): val(0) {}

    void update(const T &src) {;}

    bool ok(const T &x) const {
        return get_equilibrium(x);
    }
};

inline bool get_equilibrium(std::pair<String,answer_motoh> x){
    return x.second.openings == x.second.closings;
}

//Filter function for RNAmotiFold alignments, makes two hashmap searches and compares the outputs. This ensures only motif matching regions are checked for motifs.
template<typename alphabet, typename pos_type, typename T>
inline bool motif_match(const Basic_Sequence<alphabet,pos_type> &seq1, const Basic_Sequence<alphabet, pos_type> &seq2, T i_seq1, T j_seq1, T i_seq2, T j_seq2){
    //Initiate hashmaps
    if (!motif_ali::init) {
        fill_hashmap(gapc::Opts::getOpts()->custom_hairpins, gapc::Opts::getOpts()->replaceH, motif_ali::HairpinHashMap, motif_basic::Hairpins, motif_basic::Hairpin_lengths);
        fill_hashmap(gapc::Opts::getOpts()->custom_bulges, gapc::Opts::getOpts()->replaceB, motif_ali::BulgeHashMap, motif_basic::Bulges, motif_basic::Bulge_lengths);
        fill_hashmap(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts()->replaceI, motif_ali::InternalHashMap, motif_basic::Internals,motif_basic::Internal_lengths);
        fill_internal_hashmaps(gapc::Opts::getOpts()->custom_internals, gapc::Opts::getOpts()->replaceI, motif_ali::InternalHashMap_fronts, motif_ali::InternalHashMap_backs, motif_basic::Internals, motif_basic::Internal_lengths);
        motif_ali::init = true;
    }
    //Create Subsequence objects 
    Basic_Subsequence<alphabet, pos_type> Subseq1 {seq1, i_seq1, j_seq1};
    Basic_Subsequence<alphabet, pos_type> Subseq2 {seq2, i_seq2, j_seq2};
    std::set<char> res;
    //All motifs have this requirement actually so everything goes into this if
    if (j_seq1 < seq1.size() - 2 && j_seq2 < seq2.size() - 2 && i_seq1 > 2 && i_seq2 > 2){ //seq.size() - 2 covers end of sequence, i_seq > 2 covers beginning, this also works even in the reverse parsing case

        //Check for the internal front if there are already any found internal backs, filtering for internal front/back distance won't work here because we need to check boundaries for every possible pair, this is done in the function though!
        if (motif_ali::back_collector.size() > 0){
            std::set<char> backside = check_for_internal_front(seq1,seq2,i_seq1,j_seq1,i_seq2,j_seq2,Subseq1,Subseq2);
            res.merge(backside);
        }

        //The sequence is parsed in reverse order, so starting from the back we need to first find the backside of any internal loop
        //We don't need to look for internal loops if there are less than 9 nucleotides in the rest of the sequence because 7 (minimal hairpin loop) + 2 (minimal internal loop) = 9
        if (i_seq1 > 8 && i_seq2 > 8){
            std::set<char> Internal_backs = check_map(Subseq1,Subseq2, motif_ali::InternalHashMap_backs);
             if (Internal_backs.size() > 0){
                //Record the found backs for later checking with the internal loop fronts (no need to include the basepairs we just collect location)
                Internal_location seq_pair {{i_seq1, j_seq1} ,{i_seq2, j_seq2}};
                motif_ali::back_collector.push_back(seq_pair);
                res.merge(Internal_backs);
            }
        }
        //Hairpin Motif recognition
        if (rnali_basepairing(seq1, i_seq1 - 1, j_seq1 + 1) && rnali_basepairing(seq2, i_seq2 - 1, j_seq2 + 1) && rnali_basepairing(seq1, i_seq1 - 2 , j_seq1 + 2) && rnali_basepairing(seq2, i_seq2 - 2, j_seq2 + 2)){ // Sowohl i+1/j+1 und i+2/j+2 müssen pairen? FIXME lonely basepair version mit nur i+1/j+1
            res.merge(check_map(Subseq1, Subseq2 , motif_ali::HairpinHashMap));
        }
        //Bulges have no filter I think ? Maybe just leave bulges out in the future. They explode the search space too much.
        res.merge(check_map(Subseq1,Subseq2,motif_ali::BulgeHashMap));
    }
    return res.size() > 0;
};

template<typename alphabet, typename pos_type, typename T>
inline std::set<char> check_for_internal_front(const Basic_Sequence<alphabet,pos_type> &seq1, const Basic_Sequence<alphabet, pos_type> &seq2, T i_seq1, T j_seq1, T i_seq2, T j_seq2, Basic_Subsequence<char,unsigned int> Subseq1,Basic_Subsequence<char, unsigned int>Subseq2){
    std::set<char> res;
    for (auto paired: motif_ali::back_collector){ //Iterate over the already found internal back parts
        if (paired.first.first - i_seq1 > 8 && paired.second.first - i_seq2 > 8){
            if (rnali_basepairing(seq1,paired.first.first-2, j_seq1+2) && rnali_basepairing(seq1,paired.first.first-1,j_seq1+1) && rnali_basepairing(seq1,paired.first.second +1,i_seq1 -1) && rnali_basepairing(seq1, paired.first.second +2, i_seq1 -2)){
                if (rnali_basepairing(seq2,paired.second.first-2, j_seq2+2) && rnali_basepairing(seq2,paired.second.first-1,j_seq2+1) && rnali_basepairing(seq2,paired.second.second+1,i_seq2-1) && rnali_basepairing(seq2, paired.second.second +2, i_seq2 -2)){
                    Basic_Subsequence<alphabet, pos_type> back_first {seq1,paired.first.first,paired.first.second};
                        Basic_Subsequence<alphabet, pos_type> back_second {seq2,paired.second.first,paired.second.second};
                            if (auto search_first = motif_ali::InternalHashMap.get_motif_set(Subseq1,back_first); search_first != motif_ali::InternalHashMap.dupe_end()){
                                if (auto search_second = motif_ali::InternalHashMap.get_motif_set(Subseq2,back_first); search_second != motif_ali::InternalHashMap.dupe_end()){
                                    std::set_intersection(search_first->second.begin(),search_first->second.end(),search_second->second.begin(),search_second->second.end(),std::inserter(res,res.begin()));
                                }
                            }
                }
            }
        }
    }
    return res;
}
