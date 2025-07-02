#pragma once
#include <ostream>
#include <vector>
#include "subsequence.hh"


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
    bool empty_;

    answer_motoh() : score(0), empty_(false){}
    answer_motoh(int init_score) : score(init_score), empty_(false) {} //Called on addition/subtraction, seems suboptimal to make new obj
                                                                      //seems subopt to me but adhering to the FoldGrammars Code style.

    answer_motoh(const int init_score, seq_vector fronts_first, seq_vector fronts_second):
    score(init_score),first_track_seqs(fronts_first),second_track_seqs(fronts_second), empty_(false){}

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
        assert(!empty_);
        return {score + addition, this->first_track_seqs, this->second_track_seqs};
    }

    answer_motoh operator-(int subtract) const {
        assert(!empty_);
        return {score - subtract,this->first_track_seqs, this->second_track_seqs};
    }

    void extra_score(){
        //this is a dummy, fill out later
    }
};

inline std::ostream &operator<<(std::ostream &ostream, const answer_motoh& ans) {
    ostream << ans.score;
    return ostream;
}

inline void empty(answer_motoh &make_empty) {
    make_empty.empty_ = true;
}

inline bool isEmpty(const answer_motoh &make_empty) {
    return make_empty.empty_;
}