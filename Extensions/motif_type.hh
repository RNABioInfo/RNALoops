#pragma once
#include <ostream>
#include <utility>
#include <vector>
#include "subsequence.hh"
#include "motif.hh"
#include "motif_ali.hh"
#include "shape.hh"

#pragma once
struct answer_motoh{
#ifdef CHECKPOINTING_INTEGRATED
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & score;
        ar & motifs;
        ar & empty_
    }
#endif
    using seq_vector = std::vector<Basic_Subsequence<char, unsigned int>>;
    int score;
    seq_vector first_track_seqs;
    seq_vector second_track_seqs;
    std::vector<std::pair<Basic_Subsequence<char,unsigned int>,std::set<char>>> founds;
    bool empty_;

    answer_motoh() : score(0), empty_(false){}
    answer_motoh(int init_score) : score(init_score), empty_(false) {} //Called on addition/subtraction, seems suboptimal to make new obj
                                                                      //seems subopt to me but adhering to the FoldGrammars Code style.

    answer_motoh(const int init_score, const std::vector<std::pair<Basic_Subsequence<char,unsigned int>,std::set<char>>>&motifs,seq_vector   first_track, seq_vector second_track)  : score(init_score),first_track_seqs(std::move(first_track)),second_track_seqs(std::move(second_track)),founds(motifs),empty_(false){}

    bool operator>(const answer_motoh &other ) const {
        return score > other.score;
    }

    bool operator<(const answer_motoh &other ) const {
        return score < other.score;
    }

    bool operator==(const answer_motoh& other) const {
        return score == other.score;
    }

    bool operator<=(const answer_motoh& other) const {
        return score <= other.score;
    }

    template<typename T>
    bool operator>(const T &other) const {
        return score > other;
    }

    template<typename T>
    bool operator<(const T &other) const {
        return score < other;
    }

    template<typename T>
    bool operator==(const T &other) const {
        return score == other;
    }

    answer_motoh operator+(int addition) const {
        assert(!this.empty_);
        return {score + addition,this->founds,this->first_track_seqs,this->second_track_seqs};
    }

    answer_motoh operator-(int subtract) const {
        assert(!this.empty_);
        return {score - subtract,this->founds,this->first_track_seqs,this->second_track_seqs};
    }

    void extra_score(){
        
    }
};

inline std::ostream &operator<<(std::ostream &ostream, const answer_motoh& ans) {
    ostream << ans.score;
    return ostream;
}

inline void append(answer_motoh& bruh, const  std::set<char>& set){
    bruh.founds.emplace_back(set);
}

inline void empty(answer_motoh &make_empty) {
    make_empty.empty_ = true;
}

inline bool isEmpty(const answer_motoh &make_empty) {
    return make_empty.empty_;
}