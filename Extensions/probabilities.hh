#ifndef PROBABILITIES_HH
#define PROBABILITIES_HH

#include <utility>
#include "rnaoptions_defaults.hh"
#include "typesRNAfolding.hh"

// needed for (shape * (mfe * pf) * pretty) --kbacktrack ing since
// the 3rd component needs to be ignored (its synoptic)
template <typename SHAPE, typename MFE, typename PFUNC>
inline bool operator==(const std::pair<SHAPE, std::pair<MFE, PFUNC> > &a,
                       const std::pair<SHAPE, std::pair<MFE, PFUNC> > &b) {
  return a.first == b.first && a.second.first == b.second.first;
}

template <typename SHAPE, typename PFUNC>
inline bool operator==(
    const std::pair<SHAPE, std::pair<answer_macrostate_mfe, PFUNC> > &a,
    const std::pair<SHAPE, std::pair<answer_macrostate_mfe, PFUNC> > &b) {
  return a.first == b.first && a.second.first.energy == b.second.first.energy;
}

// needed for (mfe * pf) * pretty --backtrack
template <typename PFUNC>
inline bool operator==(const std::pair<answer_macrostate_mfe, PFUNC> &a,
                       const std::pair<answer_macrostate_mfe, PFUNC> &b) {
  return a.first == b.first;
}

template <typename SHAPE, typename MFE, typename PFUNC, typename I>
inline MFE get_pk_fn(
    const Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &candidates) {
  Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &hash =
      const_cast<Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I> &>(
          candidates);
  typename Hash::Ref<std::pair<SHAPE, std::pair<MFE, PFUNC> >, I>::iterator it =
      hash.ref().begin();
  MFE mfe;
  if (it == hash.ref().end()) {
    empty(mfe);
  } else {
    mfe = (*it).second.first;
    ++it;
    for (; it != hash.ref().end(); ++it) {
      MFE cand = (*it).second.first;
      if (getIntScore(cand) < getIntScore(mfe)) {
        mfe = cand;
      }
    }
  }
  return mfe;
}

template <typename T>
struct filterLowProbShapes {
#ifdef CHECKPOINTING_INTEGRATED
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &sum;
  }
#endif
  double sum;
  filterLowProbShapes() : sum(0) {}
  void update(const T &src) { sum += getPfuncValue(src); }
  bool ok(const T &x) const {
    double thresh = lowProbabilityFilter() * sum;
    return getPfuncValue(x) > thresh;
  }
};

template <typename SHAPE>
inline double getPfuncValue(std::pair<SHAPE, answer_macrostate_pfunc> x) {
  return sum_elems(x.second.pf);
}

template <typename SHAPE, typename PFUNC>
inline double getPfuncValue(std::pair<SHAPE, PFUNC> x) {
  return x.second;
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(std::pair<SHAPE, std::pair<MFE, double> > x) {
  return x.second.second;
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(
    std::pair<SHAPE, std::pair<MFE, answer_pknot_pfunc> > x) {
  return x.second.second.pfunc;
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(
    std::pair<SHAPE, std::pair<MFE, answer_macrostate_pfunc> > x) {
  return sum_elems(x.second.second.pf);
}

template <typename SHAPE, typename MFE>
inline double getPfuncValue(
    std::pair<SHAPE, std::pair<MFE, answer_ali_pfunc_macrostate> > x) {
  return sum_elems(x.second.second.pf);
}

#endif  // PROBABILITIES_HH
