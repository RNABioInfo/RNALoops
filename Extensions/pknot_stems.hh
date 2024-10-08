
// A dynamic programming evaluator generated by GAP-C.
//
//   GAP-C version:
//     2012-12-14 14:52 +0100 45615e44856924218a45ed3faa58059f1320c98b tip wd
//     modified
//
//   GAP-C call:
//     gapc stacklen.gap --tab-all -o stacklen_window.cc --window-mode
//
//

#ifndef pknot_stems_hh
#define pknot_stems_hh

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "rtlib/adp.hh"

typedef Basic_Subsequence<char, unsigned> TUSubsequence;

#include <rtlib/subopt.hh>

#include "rna.hh"

#define KNOT_ANSWER_TYPE answer_pknot_mfe

#ifdef WINDOW_MODE
#ifdef WITH_RNAOPTIONS
#include "rnaoptions.hh"
#endif
#include "rnaoptions_defaults.hh"
class stacklen_window {
 public:
  Basic_Sequence<char> t_0_seq;
  unsigned int t_0_left_most;
  unsigned int t_0_right_most;
  unsigned wsize;
  unsigned winc;

  std::pair<int, int> Bint_firstG_int_secondG_E_zero;

  class stack_table_t {
   private:
    unsigned wsize;
    unsigned winc;
    unsigned int t_0_left_most;
    unsigned int t_0_right_most;
    std::vector<std::pair<int, int> > array;
    std::vector<bool> tabulated;
    unsigned int t_0_n;
    std::pair<int, int> zero;
    unsigned int size() { return ((((wsize * (wsize + 1)) / 2) + wsize) + 1); }

   public:
    stack_table_t() { empty(zero); }

    void init(unsigned int t_0_n_, unsigned wsize_, unsigned winc_,
              const std::string &tname) {
      t_0_n = t_0_n_;
      t_0_left_most = 0;
      t_0_right_most = t_0_n;
      wsize = wsize_;
      winc = winc_;
      t_0_right_most = wsize;
      unsigned int newsize = size();
      array.resize(newsize);
      tabulated.clear();
      tabulated.resize(newsize);
    }
    void un_tabulate(unsigned int t_0_i, unsigned int t_0_j) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      t_0_i = (t_0_i % (wsize + 1));
      t_0_j = (t_0_j % (wsize + 1));
      if ((t_0_i > t_0_j)) {
        swap(t_0_i, t_0_j);
      }

      tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = false;
    }

    void window_increment() {
      unsigned inc = winc;
      if (t_0_left_most + winc > t_0_n) {
        inc = std::min(t_0_n - t_0_left_most, winc);
        assert(inc);
      }
      for (unsigned i = t_0_left_most; i < t_0_left_most + inc; ++i)
        for (unsigned j = i; j <= t_0_right_most; ++j) {
          un_tabulate(i, j);
        }
      t_0_left_most += inc;
      t_0_right_most = std::min(t_0_right_most + inc, t_0_n);
    }

    bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      t_0_i = (t_0_i % (wsize + 1));
      t_0_j = (t_0_j % (wsize + 1));
      if ((t_0_i > t_0_j)) {
        swap(t_0_i, t_0_j);
      }

      return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
    }

    void clear() { tabulated.clear(); }
    std::pair<int, int> &get(unsigned int t_0_i, unsigned int t_0_j) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      t_0_i = (t_0_i % (wsize + 1));
      t_0_j = (t_0_j % (wsize + 1));
      if ((t_0_i > t_0_j)) {
        swap(t_0_i, t_0_j);
      }

      assert(tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
      assert(((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
      return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
    }

    void set(unsigned int t_0_i, unsigned int t_0_j,
             const std::pair<int, int> &e) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      assert(!is_tabulated(t_0_i, t_0_j));
      t_0_i = (t_0_i % (wsize + 1));
      t_0_j = (t_0_j % (wsize + 1));
      if ((t_0_i > t_0_j)) {
        swap(t_0_i, t_0_j);
      }

      assert(((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
      array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
      tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
    }
  };
  stack_table_t stack_table;

  template <typename C>
  void init(const Basic_Sequence<C> &seq, unsigned int window_size,
            unsigned int window_increment) {
    t_0_seq = seq;
    stack_table.init(t_0_seq.size(), window_size, window_increment,
                     "stack_table");
    empty(Bint_firstG_int_secondG_E_zero);

    t_0_left_most = 0;
    t_0_right_most = t_0_seq.size();
    t_0_right_most = window_size;
    wsize = window_size;
    winc = window_increment;
  }

  void window_increment() {
    stack_table.window_increment();
    t_0_left_most += winc;
    t_0_right_most = std::min(t_0_seq.size(), t_0_left_most + wsize);
  }

 public:
  std::pair<int, int> &nt_stack(unsigned int t_0_i, unsigned int t_0_j) {
    if (stack_table.is_tabulated(t_0_i, t_0_j)) {
      return stack_table.get(t_0_i, t_0_j);
    }

    List_Ref<std::pair<int, int> > answers;
    empty(answers);
    empty(answers);
    std::pair<int, int> ret_0;
    if (((t_0_j - t_0_i) >= 2)) {
      if (basepairing(t_0_seq, t_0_i, t_0_j)) {
        TUSubsequence ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j);
        TUSubsequence a_2 = ret_3;
        if (is_not_empty(a_2)) {
          TUSubsequence ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1));
          TUSubsequence a_0 = ret_1;
          if (is_not_empty(a_0)) {
            std::pair<int, int> ret_2 = nt_stack((t_0_i + 1), (t_0_j - 1));
            std::pair<int, int> a_1 = ret_2;
            if (is_not_empty(a_1)) {
              ret_0 = sr(a_0, a_1, a_2);
            } else {
              empty(ret_0);
            }

            erase(a_1);
          } else {
            empty(ret_0);
          }

          erase(a_0);
        } else {
          empty(ret_0);
        }

        erase(a_2);
      } else {
        empty(ret_0);
        empty(ret_0);
      }
    } else {
      empty(ret_0);
    }

    if (is_not_empty(ret_0)) {
      push_back_max_other(answers, ret_0);
    }

    std::pair<int, int> ret_4;
    if (((t_0_j - t_0_i) >= 0)) {
      TUSubsequence ret_5 = REGION0(t_0_seq, t_0_i, t_0_j);
      TUSubsequence a_3 = ret_5;
      if (is_not_empty(a_3)) {
        ret_4 = end(a_3);
      } else {
        empty(ret_4);
      }

      erase(a_3);
    } else {
      empty(ret_4);
    }

    if (is_not_empty(ret_4)) {
      push_back_max_other(answers, ret_4);
    }

    std::pair<int, int> eval = h(answers);
    erase(answers);
    stack_table.set(t_0_i, t_0_j, eval);
    return stack_table.get(t_0_i, t_0_j);
  }

 private:
  std::pair<int, int> end(const TUSubsequence &p_e) {
    TUSubsequence l_0 = p_e;
    TUSubsequence r_0 = p_e;
    int ret_left = end_l(l_0);
    int ret_right = end_r(r_0);
    std::pair<int, int> ret;
    ret.first = ret_left;
    ret.second = ret_right;
    return ret;
  }

  std::pair<int, int> h(List_Ref<std::pair<int, int> > i) {
    std::pair<List<std::pair<int, int> >::iterator,
              List<std::pair<int, int> >::iterator>
        range = get_range(i);
    return h(range);
  }

  template <typename Iterator>
  std::pair<int, int> h(std::pair<Iterator, Iterator> i) {
    std::pair<int, int> answers;
    empty(answers);
    std::pair<
        Proxy::Iterator<Iterator, select1st<typename Iterator::value_type> >,
        Proxy::Iterator<Iterator, select1st<typename Iterator::value_type> > >
        left = splice_left(i);
    int left_answers = h_l(left);
    if (isEmpty(left_answers)) {
      std::pair<int, int> temp;
      empty(temp);
      erase(left_answers);
      return temp;
    }

    int elem = left_answers;
    List_Ref<int> right_candidates;
    empty(right_candidates);
    empty(right_candidates);
    for (Iterator tupel = i.first; tupel != i.second; ++tupel) {
      if (((*tupel).first == elem)) {
        push_back(right_candidates, (*tupel).second);
      }
    }
    int right_answers = h_r(right_candidates);
    int right_elem;
    right_elem = right_answers;
    std::pair<int, int> temp_elem;
    temp_elem.first = elem;
    temp_elem.second = right_elem;
    answers = temp_elem;
    return answers;
  }

  std::pair<int, int> sr(const TUSubsequence &p_l,
                         const std::pair<int, int> &p_x,
                         const TUSubsequence &p_r) {
    TUSubsequence l_0 = p_l;
    int l_1 = p_x.first;
    TUSubsequence l_2 = p_r;
    TUSubsequence r_0 = p_l;
    int r_1 = p_x.second;
    TUSubsequence r_2 = p_r;
    int ret_left = sr_l(l_0, l_1, l_2);
    int ret_right = sr_r(r_0, r_1, r_2);
    std::pair<int, int> ret;
    ret.first = ret_left;
    ret.second = ret_right;
    return ret;
  }

  int end_l(const TUSubsequence &e) { return 0; }

  int h_l(List_Ref<int> i) {
    std::pair<List<int>::iterator, List<int>::iterator> range = get_range(i);
    return h_l(range);
  }

  template <typename Iterator>
  int h_l(std::pair<Iterator, Iterator> i) {
    return maximum(i);
  }

  int sr_l(const TUSubsequence &l, int x, const TUSubsequence &r) {
    return (1 + x);
  }

  int end_r(const TUSubsequence &e) { return 0; }

  int h_r(List_Ref<int> i) {
    std::pair<List<int>::iterator, List<int>::iterator> range = get_range(i);
    return h_r(range);
  }

  template <typename Iterator>
  int h_r(std::pair<Iterator, Iterator> i) {
    return minimum(i);
  }

  int sr_r(const TUSubsequence &l, int x, const TUSubsequence &r) {
    return (x + sr_energy(l, r));
  }

 public:
  void cyk() {}

 public:
  std::pair<int, int> &run() { return nt_stack(t_0_left_most, t_0_right_most); }
  void print_stats(std::ostream &o) {
#ifdef STATS
    o << "\n\nN = " << seq.size() << '\n';
    stack_table.print_stats(o, "stack_table");
#endif
  }

  template <typename Value>
  void print_result(std::ostream &out, Value &res) {
    if (isEmpty(res))
      out << "[]\n";
    else
      out << res << '\n';
  }
  template <typename Value>
  void print_backtrack(std::ostream &out, Value &value)

  {}
  void print_subopt(std::ostream &out, int delta = 0) {}
};

template <typename C, typename U>
inline std::pair<int, unsigned> stacklen(const Basic_Sequence<C> &seq, U a,
                                         U b) {
  static stacklen_window pkStems;
  static std::map<U, std::map<U, bool> > isComputed;

  if (isComputed[a][b]) {
    return pkStems.nt_stack(a, b);
  } else {
    pkStems.init(seq, getWindowSize(), getWindowIncrement());
    pkStems.run();
    isComputed[a][b] = true;
    return pkStems.nt_stack(a, b);
  }
}

template <typename A, typename B>
inline B &energy(std::pair<A, B> &p) {
  return p.second;
}

template <typename A, typename B>
inline A &length(std::pair<A, B> &p) {
  return p.first;
}

template <typename A, typename B>
inline const B &energy(const std::pair<A, B> &p) {
  return p.second;
}

template <typename A, typename B>
inline const A &length(const std::pair<A, B> &p) {
  return p.first;
}

#else
class stacklen_normal {
 public:
  Basic_Sequence<char> t_0_seq;
  unsigned int t_0_left_most;
  unsigned int t_0_right_most;

  std::pair<int, int> Bint_firstG_int_secondG_E_zero;

  class stack_table_t {
   private:
    unsigned int t_0_left_most;
    unsigned int t_0_right_most;
    std::vector<std::pair<int, int> > array;
    std::vector<bool> tabulated;
    unsigned int t_0_n;
    std::pair<int, int> zero;
    unsigned int size() {
      return (1 * ((((t_0_n * (t_0_n + 1)) / 2) + t_0_n) + 1));
    }

   public:
    stack_table_t() { empty(zero); }

    void init(unsigned int t_0_n_, const std::string &tname) {
      t_0_n = t_0_n_;
      t_0_left_most = 0;
      t_0_right_most = t_0_n;
      unsigned int newsize = size();
      array.resize(newsize);
      tabulated.clear();
      tabulated.resize(newsize);
    }
    bool is_tabulated(unsigned int t_0_i, unsigned int t_0_j) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      return tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
    }

    void clear() { tabulated.clear(); }
    std::pair<int, int> &get(unsigned int t_0_i, unsigned int t_0_j) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      assert(tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))]);
      assert(((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
      return array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))];
    }

    void set(unsigned int t_0_i, unsigned int t_0_j,
             const std::pair<int, int> &e) {
      assert((t_0_i <= t_0_j));
      assert((t_0_j <= t_0_n));
      assert(!is_tabulated(t_0_i, t_0_j));
      assert(((0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i))) < size()));
      array[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = e;
      tabulated[(0 + (1 * (((t_0_j * (t_0_j + 1)) / 2) + t_0_i)))] = true;
    }
  };
  stack_table_t stack_table;

  template <typename C>
  void init(const Basic_Sequence<C> &seq) {
    t_0_seq = seq;
    t_0_seq.n = seq.n;

    stack_table.init(t_0_seq.size(), "stack_table");
    empty(Bint_firstG_int_secondG_E_zero);
    t_0_left_most = 0;
    t_0_right_most = t_0_seq.size();
  }

  // void init(const gapc::Opts &opts)
  //{
  // const std::vector<std::pair<const char *, unsigned> > &inp = opts.inputs;
  // if(inp.size() != 1)
  //   throw gapc::OptException("Number of input sequences does not match.");
  //
  //   t_0_seq.copy(inp[0].first, inp[0].second);
  // char_to_rna(t_0_seq);
  //   stack_table.init( t_0_seq.size(), "stack_table");
  // empty(Bint_firstG_int_secondG_E_zero);
  //
  // t_0_left_most = 0;
  // t_0_right_most = t_0_seq.size();
  // }

 public:
  std::pair<int, int> &nt_stack(unsigned int t_0_i, unsigned int t_0_j) {
    if (stack_table.is_tabulated(t_0_i, t_0_j)) {
      return stack_table.get(t_0_i, t_0_j);
    }

    List_Ref<std::pair<int, int> > answers;
    empty(answers);
    empty(answers);
    std::pair<int, int> ret_0;
    if (((t_0_j - t_0_i) >= 2)) {
      if (basepairing(t_0_seq, t_0_i, t_0_j)) {
        TUSubsequence ret_3 = BASE(t_0_seq, (t_0_j - 1), t_0_j);
        TUSubsequence a_2 = ret_3;
        if (is_not_empty(a_2)) {
          TUSubsequence ret_1 = BASE(t_0_seq, t_0_i, (t_0_i + 1));
          TUSubsequence a_0 = ret_1;
          if (is_not_empty(a_0)) {
            std::pair<int, int> ret_2 = nt_stack((t_0_i + 1), (t_0_j - 1));
            std::pair<int, int> a_1 = ret_2;
            if (is_not_empty(a_1)) {
              ret_0 = sr(a_0, a_1, a_2);
            } else {
              empty(ret_0);
            }

            erase(a_1);
          } else {
            empty(ret_0);
          }

          erase(a_0);
        } else {
          empty(ret_0);
        }

        erase(a_2);
      } else {
        empty(ret_0);
        empty(ret_0);
      }

    } else {
      empty(ret_0);
    }

    if (is_not_empty(ret_0)) {
      push_back_max_other(answers, ret_0);
    }

    std::pair<int, int> ret_4;
    if (((t_0_j - t_0_i) >= 0)) {
      TUSubsequence ret_5 = REGION0(t_0_seq, t_0_i, t_0_j);
      TUSubsequence a_3 = ret_5;
      if (is_not_empty(a_3)) {
        ret_4 = end(a_3);
      } else {
        empty(ret_4);
      }

      erase(a_3);
    } else {
      empty(ret_4);
    }

    if (is_not_empty(ret_4)) {
      push_back_max_other(answers, ret_4);
    }

    std::pair<int, int> eval = h(answers);
    erase(answers);
    stack_table.set(t_0_i, t_0_j, eval);
    return stack_table.get(t_0_i, t_0_j);
  }

 private:
  std::pair<int, int> end(const TUSubsequence &p_e) {
    TUSubsequence l_0 = p_e;
    TUSubsequence r_0 = p_e;
    int ret_left = end_l(l_0);
    int ret_right = end_r(r_0);
    std::pair<int, int> ret;
    ret.first = ret_left;
    ret.second = ret_right;
    return ret;
  }

  std::pair<int, int> h(List_Ref<std::pair<int, int> > i) {
    std::pair<List<std::pair<int, int> >::iterator,
              List<std::pair<int, int> >::iterator>
        range = get_range(i);
    return h(range);
  }

  template <typename Iterator>
  std::pair<int, int> h(std::pair<Iterator, Iterator> i) {
    std::pair<int, int> answers;
    empty(answers);
    std::pair<
        Proxy::Iterator<Iterator, select1st<typename Iterator::value_type> >,
        Proxy::Iterator<Iterator, select1st<typename Iterator::value_type> > >
        left = splice_left(i);
    int left_answers = h_l(left);
    if (isEmpty(left_answers)) {
      std::pair<int, int> temp;
      empty(temp);
      erase(left_answers);
      return temp;
    }

    int elem = left_answers;
    List_Ref<int> right_candidates;
    empty(right_candidates);
    empty(right_candidates);
    for (Iterator tupel = i.first; tupel != i.second; ++tupel) {
      if (((*tupel).first == elem)) {
        push_back(right_candidates, (*tupel).second);
      }
    }
    int right_answers = h_r(right_candidates);
    int right_elem;
    right_elem = right_answers;
    std::pair<int, int> temp_elem;
    temp_elem.first = elem;
    temp_elem.second = right_elem;
    answers = temp_elem;
    return answers;
  }

  std::pair<int, int> sr(const TUSubsequence &p_l,
                         const std::pair<int, int> &p_x,
                         const TUSubsequence &p_r) {
    TUSubsequence l_0 = p_l;
    int l_1 = p_x.first;
    TUSubsequence l_2 = p_r;
    TUSubsequence r_0 = p_l;
    int r_1 = p_x.second;
    TUSubsequence r_2 = p_r;
    int ret_left = sr_l(l_0, l_1, l_2);
    int ret_right = sr_r(r_0, r_1, r_2);
    std::pair<int, int> ret;
    ret.first = ret_left;
    ret.second = ret_right;
    return ret;
  }

  int end_l(const TUSubsequence &e) { return 0; }
  int h_l(List_Ref<int> i) {
    std::pair<List<int>::iterator, List<int>::iterator> range = get_range(i);
    return h_l(range);
  }

  template <typename Iterator>
  int h_l(std::pair<Iterator, Iterator> i) {
    return maximum(i);
  }

  int sr_l(const TUSubsequence &l, int x, const TUSubsequence &r) {
    return (1 + x);
  }

  int end_r(const TUSubsequence &e) { return 0; }
  int h_r(List_Ref<int> i) {
    std::pair<List<int>::iterator, List<int>::iterator> range = get_range(i);
    return h_r(range);
  }
  template <typename Iterator>
  int h_r(std::pair<Iterator, Iterator> i) {
    return minimum(i);
  }
  int sr_r(const TUSubsequence &l, int x, const TUSubsequence &r) {
    return (x + sr_energy(l, r));
  }

 public:
  void cyk();

 public:
  std::pair<int, int> &run() { return nt_stack(t_0_left_most, t_0_right_most); }
  void print_stats(std::ostream &o) {
#ifdef STATS
    o << "\n\nN = " << seq.size() << '\n';
    stack_table.print_stats(o, "stack_table");
#endif
  }

  template <typename Value>
  void print_result(std::ostream &out, Value &res) {
    if (isEmpty(res))
      out << "[]\n";
    else
      out << res << '\n';
  }
  template <typename Value>
  void print_backtrack(std::ostream &out, Value &value)

  {}
  void print_subopt(std::ostream &out, int delta = 0) {}
};

template <typename C, typename U>
inline std::pair<int, int> stacklen(const Basic_Sequence<C> &seq, U a, U b) {
  static stacklen_normal pkStems;
  static bool compute = true;

  if (!compute) {
    return pkStems.nt_stack(a, b);
  } else {
    pkStems.init(seq);
    pkStems.run();
    compute = false;
    return pkStems.nt_stack(a, b);
  }
}

template <typename A, typename B>
inline A &energy(std::pair<A, B> &p) {
  return p.second;
}

template <typename A, typename B>
inline B &length(std::pair<A, B> &p) {
  return p.first;
}

template <typename A, typename B>
inline const A &energy(const std::pair<A, B> &p) {
  return p.second;
}

template <typename A, typename B>
inline const B &length(const std::pair<A, B> &p) {
  return p.first;
}

#endif

#endif
