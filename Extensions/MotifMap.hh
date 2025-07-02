#pragma once
#include <algorithm>
#include<iostream>
#include <numeric>
#include <optional>
#include<set>
#include <string>
#include<unordered_map>
#include<string_view>
#include<utility>
#include "ali_t.hh"
#include "sequence.hh"
#include "subsequence.hh"


struct std_set_hash {
    template<typename T1>
    size_t operator () (const std::set<T1> &input_set) const {
        std::size_t seed = input_set.size();
        for (const auto& iter: input_set) {
            seed ^= iter + 0x9e3779b9  + (seed << 6) + (seed >> 2); //NOLINT
        }
        return seed;
    }
};

struct Motif {std::string seq; char abb {'X'};};

enum class direction_type: std::uint8_t {forward,reverse,both};

class MotifMap{
    
    public:
    using DupeHashMap = std::unordered_map<Basic_Sequence<char, unsigned int>, std::set<char>, Hash_ali_array>; //Type for the Dupe Map recording all possible motifs for each sequence
    using MotifHashMap = std::unordered_map<Basic_Sequence<char, unsigned int>, char, Hash_ali_array>; //Type for the central Motif Map that contains motifs with their abbreviations
    using TranslationHashMap = std::unordered_map<std::set<char>, char,std_set_hash>; //Type for the translation map that connects 
    const direction_type direction;
    std::string_view input;
    std::set<char> abbreviations;

    //Constructor used for direction only input, used mostly for ali cases probably
    MotifMap(direction_type dir): direction(dir){} 

    //Most used constructor to create a MotifHash Map from an input string view and a direction for the motifs
    MotifMap(const std::string_view input_string, direction_type dir): direction(dir),input(input_string){ //Constructor definition
        fill_Dupes(input, direction);
        fill_Translations();
        fill_Motifs();
    };

    template<typename alphabet, typename pos_type>
    const DupeHashMap::iterator get_motif_set(const Basic_Subsequence<alphabet, pos_type> &input_subsequence){
        Basic_Sequence Motif {&input_subsequence.front(), input_subsequence.size()};
        return Dupes.find(Motif);
    };

    template<typename alphabet, typename pos_type>
    const MotifHashMap::iterator get_motif(const Basic_Subsequence<alphabet, pos_type> &input_subsequence){
        Basic_Sequence Motif {&input_subsequence.front(),input_subsequence.size()};
        return Motifs.find(Motif);
    };

    template<typename alphabet, typename pos_type>
    const DupeHashMap::iterator get_motif_set(const Basic_Sequence<alphabet, pos_type> &input_sequence){
        return Dupes.find(input_sequence);
    };

    template<typename alphabet, typename pos_type>
    const MotifHashMap::iterator get_motif(const Basic_Sequence<alphabet,pos_type> &input_sequence){
        return Motifs.find(input_sequence);
    };

    template<typename alphabet, typename pos_type>
    const MotifHashMap::iterator get_motif(const Basic_Subsequence<alphabet, pos_type> &internal_subsequence1, const Basic_Subsequence<alphabet, pos_type> &internal_subsequence2){
        Basic_Sequence Motif1 {&internal_subsequence1.front(),internal_subsequence1.size()};
        Basic_Sequence Motif2 {&internal_subsequence2.front(),internal_subsequence2.size()};
        Motif1.concat(Motif2.seq,Motif2.size());
        return Motifs.find(Motif1);
    };

    template<typename alphabet, typename pos_type>
    const DupeHashMap::iterator get_motif_set(const Basic_Subsequence<alphabet, pos_type> &internal_subsequence1, const Basic_Subsequence<alphabet, pos_type> &internal_subsequence2){
        Basic_Sequence Motif1 {&internal_subsequence1.front(),internal_subsequence1.size()};
        Basic_Sequence Motif2 {&internal_subsequence2.front(),internal_subsequence2.size()};
        Motif1.concat(Motif2.seq,Motif2.size());
        return Dupes.find(Motif1);
    };

    template<typename alphabet, typename pos_type>
    const std::set<char> get_motif_set_ali(const Basic_Subsequence<alphabet, pos_type> &align_seq1, const Basic_Subsequence<alphabet, pos_type> &align_seq2){
        Basic_Sequence Motif1 {align_seq1.front(),align_seq1.size()};
        Basic_Sequence Motif2 {align_seq2.front(),align_seq2.size()};
        std::set<char> returnset {};
        if (auto search1 = Dupes.find(Motif1); search1 != Dupes.end()){
            returnset.merge(search1->second);
        }
        if (auto search2 = Dupes.find(Motif2); search2 != Dupes.end()){
            returnset.merge(search2->second);
        }
        return returnset;
    }

    MotifHashMap::iterator begin() {return Motifs.begin();};
    MotifHashMap::iterator end() {return Motifs.end();};

    DupeHashMap::iterator dupe_begin() {return Dupes.begin();};
    DupeHashMap::iterator dupe_end() {return Dupes.end();};

    DupeHashMap::const_iterator dupe_begin() const {return Dupes.begin();};
    DupeHashMap::const_iterator dupe_end() const {return Dupes.end();};

    MotifHashMap::const_iterator begin() const {return Motifs.begin();};
    MotifHashMap::const_iterator end() const {return Motifs.end();};

    void add_motifs(std::string_view new_input){
        switch (direction) {
            case direction_type::forward:
                add_dupe_entries(new_input, false);
                break;
            case direction_type::reverse:
                add_dupe_entries(new_input, true);
                break;
            case direction_type::both:
                add_dupe_entries(new_input, false);
                add_dupe_entries(new_input, true);
                break;
        }
    };

    void update_motif_hashmap(){
        fill_Translations();
        fill_Motifs();
    }

    void print_duplicates(){
        std::vector<std::string> strings;
        for (const auto &keyvalue: Translations){
            std::vector<std::string> substring_vector;
            if (keyvalue.first.size() > 1){
            for (char mot: keyvalue.first){
                std::string substring(1,mot);
                substring_vector.emplace_back(substring);
            }
            std::string sub = std::accumulate(std::next(substring_vector.begin()), substring_vector.end(),substring_vector[0],[](const std::string& Existing_string, const std::string& Added_string){ return Existing_string +"/" + Added_string;});
            sub += " -> ";
            sub += keyvalue.second;
            strings.emplace_back(sub);
        }}
        if (strings.size() > 0) {
            std::string dupe_string = std::accumulate(std::next(strings.begin()), strings.end(),strings[0],[](const std::string& Existing_string, const std::string& Added_string){ return Existing_string +", " + Added_string;});
            std::cerr << "Attention, the same sequence appears for different motifs. Ambiguity cases: " << dupe_string << "\n";
        }
    };

    private:
    DupeHashMap Dupes;
    TranslationHashMap Translations;
    MotifHashMap Motifs;
    

    //Initializer functions
    void fill_Dupes(const std::string_view input, direction_type direction){
        switch (direction) {
            case direction_type::forward:
                add_dupe_entries(input, false);
                break;
            case direction_type::reverse:
                add_dupe_entries(input, true);
                break;
            case direction_type::both:
                add_dupe_entries( input, false);
                add_dupe_entries( input, true);
                break;
        }
    };

    void fill_Translations(){
        std::set<std::set<char>> dupe_sets;
        std::set<char> used_characters;
        for (const auto& keyvalue: Dupes) {
            dupe_sets.insert(keyvalue.second);
        }
        for (const std::set<char>& set: dupe_sets){
            bool inserted = false;
            for (const char option: set){
                if (set.size() == 1 ){
                    used_characters.insert(option);
                    Translations[set] = option;
                    inserted = true;
                    break;
                }
                if (used_characters.find(std::tolower(option,std::locale())) == used_characters.end()){
                    used_characters.insert(std::tolower(option,std::locale()));
                    Translations[set] = std::tolower(option,std::locale());
                    inserted = true;
                    break;
                }
                continue;
            if (!inserted){
                throw std::runtime_error("One of the sets can not be used because all characters in it are used for other sets");
            }
            }
            //This is where some sort of logic would be nice to switch around the keys to see if any other combination is possible
        }
    };

    void fill_Motifs(){
        for (const auto& [key,value]: Dupes){
            Motifs[key] = Translations.at(value);
        }
    };

    //General purpose in class functions
    //Function to add entries into DupeHashMap 
    void add_dupe_entries(const std::string_view input_string_v, bool dir_bool){
        size_t pos_b = 0;
        size_t pos_e = 0;
        Basic_Sequence<char,unsigned int> basic_motif;
        while ((pos_e = input_string_v.find('\n',pos_b)) != std::string::npos) {
            split_parse_line(pos_b, pos_e, input_string_v, dir_bool);
            pos_b = pos_e + 1 ;
        }
        if (pos_b <= input_string_v.size()){ //This catches the edge case of the last line not having a new line to use as a separator
            split_parse_line(pos_b, input_string_v.size()-pos_b, input_string_v,  dir_bool);
        }
    };

    void split_parse_line(const size_t begin, const size_t end, const std::string_view view, bool dir_bool){
        std::string line;
        line = view.substr(begin,end-begin);
        if (const auto motif = parse_motif(line, dir_bool)){
            Basic_Sequence<char, unsigned int> basic_seq_motif {motif.value().seq.data(), static_cast<unsigned int>(motif.value().seq.size())};
            char_to_rnali(basic_seq_motif);
            if (auto dupe_sarch = Dupes.find(basic_seq_motif); dupe_sarch == Dupes.end()){ //Check if the sequence is already in the map
                Dupes[basic_seq_motif] = std::set<char> {motif.value().abb};    
            }
            else {
                Dupes[basic_seq_motif].insert(motif.value().abb);
            }
            abbreviations.insert(motif.value().abb);
        }
    };

    //Parsing function for input strings, skips over empty lines.
    static std::optional<Motif> parse_motif (const std::string &input, bool rev) {
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
                std::cerr << "Error more than one motif in one line\n";
                return std::nullopt;
            }
            ++elements;
        }
        if (elements > 2){
            std::cerr << "Found " << elements << ", Expected " << 2 << "\n";
            return std::nullopt;
        }
        return motif;
    };
};