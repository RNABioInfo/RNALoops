#ifndef RNAOPTIONS_DEFAULTS_HH
#define RNAOPTIONS_DEFAULTS_HH

#include <limits>
#include <utility>
#include <vector>
#include <boost/math/special_functions/fpclassify.hpp>  // isnan

#ifdef WITH_RNAOPTIONS
// use command line parameter options to define energy penalties for
// initializing pseudoknots, minimal length of kissing hairpin stems and the
// pKiss strategy
#include "rnaoptions.hh"

inline static int alignment_match() {
  return gapc::Opts::getOpts()->match_score;
}
inline static int alignment_mismatch() {
  return gapc::Opts::getOpts()->mismatch_score;
}
inline static int alignment_gap_open() {
  return gapc::Opts::getOpts()->gap_open;
}
inline static int alignment_gap_extension() {
  return gapc::Opts::getOpts()->gap_extension;
}
inline static int pkinit() {  // initialization cost for opening a new pseudoknot. Default is 900.
  return gapc::Opts::getOpts()->energyPenaltyHtype;
}
inline static int pkissinit() {  // initialization cost for opening a new
                                 // kissing hairpin. Default is 1200.
  return gapc::Opts::getOpts()->energyPenaltyKtype;
}
inline static int
minLengthKissingHairpinStems() {  // minimal length of those two stems in a KH
                                  // that form the hairpins, not the crossing
                                  // stem of the kiss. Default is 2
  return gapc::Opts::getOpts()->minimalHelixLength;
}
inline static int maxPseudoknotSize() {
  return gapc::Opts::getOpts()->maximalPseudoknotSize;
}
inline static float
lowProbabilityFilter() {  // heuristically filtering out shapes with in very low
                          // initial probability. Default is 10^-6
  return gapc::Opts::getOpts()->lowProbabilityFilter;
}
inline static int shapelevel() { return gapc::Opts::getOpts()->shapelevel; }
template <typename alphabet, typename pos_type, typename T>
inline bool selectStrategy(const Basic_Sequence<alphabet, pos_type>& seq, T i,
                           T j, const char strategy) {
  return gapc::Opts::getOpts()->strategy == strategy;
}
template <typename alphabet, typename pos_type, typename T>
inline bool allowLonelyBasepairs(const Basic_Sequence<alphabet, pos_type>& seq,
                                 T i, T j, const bool isLonelyBP) {
  return gapc::Opts::getOpts()->allowLonelyBasepairs == isLonelyBP;
}
inline int getSuboptRange(
    int mfe) {  // use command line parameter options to define the range of
                // suboptimal answers, depending on MFE.
  int range = mfe + static_cast<int>(
                        gapc::Opts::getOpts()->energydeviation_absolute * 100);
  if ((boost::math::isnan)(gapc::Opts::getOpts()->energydeviation_absolute)) {
    range = mfe *
            (100 - gapc::Opts::getOpts()->energydeviation_relative *
                       (mfe < 0 ? 1 : -1)) /
            100;
  }
  return range;
}
inline unsigned int getWindowSize() {
  return gapc::Opts::getOpts()->window_size;
}
inline unsigned int getWindowIncrement() {
  return gapc::Opts::getOpts()->window_increment;
}
inline static float getAlifold_cfactor() {
  return gapc::Opts::getOpts()->alifold_cfactor;
}
inline static float getAlifold_nfactor() {
  return gapc::Opts::getOpts()->alifold_nfactor;
}
inline static float getAlifold_minscore_basepair() {
  return gapc::Opts::getOpts()->alifold_minscore_basepair;
}
inline static const char* getDotplotFilename() {
  return gapc::Opts::getOpts()->dotPlotFilename;
}
inline static const char* getProbing_dataFilename() {
  return gapc::Opts::getOpts()->probing_dataFilename;
}
inline static float getProbing_slope() {
  return gapc::Opts::getOpts()->probing_slope;
}
inline static float getProbing_intercept() {
  return gapc::Opts::getOpts()->probing_intercept;
}
inline static const char* getProbing_modifier() {
  return gapc::Opts::getOpts()->probing_modifier;
}
inline static const char* getProbing_normalization() {
  return gapc::Opts::getOpts()->probing_normalization;
}
inline static std::vector<std::pair<const char*, unsigned> > getInputs() {
  return gapc::Opts::getOpts()->inputs;
}
inline static int getConsensusType() {
  return gapc::Opts::getOpts()->consensusType;
}
inline static bool isRiboseScoring() {
  return gapc::Opts::getOpts()->ribosum_scoring;
}
#else
// if compiled with no special options to ask for energy penalties for
// initializing pseudoknots, minimal length of kissing hairpin stems and the
// pKiss strategy.
inline static int
pkinit() {  // initialization cost for opening a new pseudoknot. Default is 900.
  return 900;
}
inline static int pkissinit() {  // initialization cost for opening a new
                                 // kissing hairpin. Default is 1200.
  return 1200;
}
inline static int
minLengthKissingHairpinStems() {  // minimal length of those two stems in a KH
                                  // that form the hairpins, not the crossing
                                  // stem of the kiss. Default is 2
  return 2;
}
inline static int maxPseudoknotSize() {
  return std::numeric_limits<int>::max();
}
inline static float lowProbabilityFilter() { return 0.000001; }
inline static int shapelevel() { return 5; }
template <typename alphabet, typename pos_type, typename T>
inline bool selectStrategy(const Basic_Sequence<alphabet, pos_type>& seq, T i,
                           T j, const char strategy) {
  return 'A' == strategy;
}
template <typename alphabet, typename pos_type, typename T>
inline bool allowLonelyBasepairs(const Basic_Sequence<alphabet, pos_type>& seq,
                                 T i, T j, const bool isLonelyBP) {
  return false == isLonelyBP;
}
inline int getSuboptRange(
    int mfe) {  // if compiled with no special options to ask for energy range,
                // use 5% of MFE as a default.
  return mfe * (100 - 5 * (mfe < 0 ? 1 : -1)) / 100;
}
inline unsigned int getWindowSize() { return 7; }
inline unsigned int getWindowIncrement() { return 1; }
inline static float getAlifold_cfactor() { return 1.0; }
inline static float getAlifold_nfactor() { return 1.0; }
inline static float getAlifold_minscore_basepair() { return -200; }
inline static const char* getDotplotFilename() { return "./dotPlot.ps"; }
inline static const char* getProbing_dataFilename() { return "\0"; }
inline static float getProbing_slope() { return 1.8 * 100; }
inline static float getProbing_intercept() { return -0.6 * 100; }
inline static const char* getProbing_modifier() { return "unknown"; }
inline static const char* getProbing_normalization() { return "centroid"; }
inline static const char* getCustom_hairpins() {return "\0";}
inline static const char* getCustom_internals() {return "\0";}
inline static const char* getCustom_bulges() {return "\0";}
inline static std::vector<std::pair<const char*, unsigned> > getInputs() {
  std::vector<std::pair<const char*, unsigned> > emptyInput;
  emptyInput.push_back(std::make_pair(const_cast<char*>(""), 0));
  return (emptyInput);
}
inline static int getConsensusType() { return 0; }
inline static bool isRiboseScoring() { return false; }
#endif

#endif  // RNAOPTIONS_DEFAULTS_HH
