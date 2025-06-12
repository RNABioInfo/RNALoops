#pragma once
#include <ostream>
#include <vector>
#include "motif.hh"
#include "motif_ali.hh"
#include "shape.hh"

#pragma once
struct answer_motoh{
#ifdef CHECKPOINTING_INTEGRATED
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        //how to serialize my struct?
    }
#endif

    int score;
    std::vector<std::set<char>> founds;
    bool empty_;

    answer_motoh() : score(0), empty_(false){}
    answer_motoh(int init_score) : score(init_score), empty_(false) {} //Called on addition/subtraction, seems suboptimal to make new obj
                                                                      //seems subopt to me but adhering to the FoldGrammars Code style.
    answer_motoh(int init_score, std::vector<std::set<char>> motifs) : score(init_score),founds(motifs),empty_(false){}

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
        return {score + addition,this->founds};
    }

    answer_motoh operator-(int subtract) const {
        assert(!empty_);
        return {score - subtract,this->founds};
    }

};

inline std::ostream &operator<<(std::ostream &o, const answer_motoh& ans) {
    o << ans.score;
    return o;
}

inline void append(answer_motoh& bruh, const  std::set<char> set){
    bruh.founds.emplace_back(set);
}

inline void empty(answer_motoh &e) {
    e.empty_ = true;
}

inline bool isEmpty(const answer_motoh &e) {
    return e.empty_;
}