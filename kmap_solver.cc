#include "kmap_solver.h"

#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <unordered_set>

#include "absl/strings/str_join.h"

namespace kmap_solver {
namespace {
// Binary increment, assuming that LSB is 0, MSB is x.Count
bool BinaryIncrement(std::vector<bool>* x) {
  std::vector<bool>& xx = *x;
  for (size_t i = 0; i < xx.size(); i++) {
    
    if (!xx[i]) {
      xx[i] = true;
      return false;
    }
    xx[i] = false;
  }
  return true;
}
template <typename T>
class SubSet {
public:
  class iterator {
  public:
    typedef iterator self_type;
    typedef T value_type;
    typedef const T& reference;
    typedef const T* pointer;
    typedef std::forward_iterator_tag iterator_category;
    typedef ptrdiff_t difference_type;

    iterator(
      typename std::vector<T>::const_iterator v,
      std::vector<bool>::const_iterator s,
      std::vector<bool>::const_iterator s_end) : v_(v), s_(s), s_end_(s_end) {}

    self_type operator++() {
      do {
        ++s_; ++v_;
      } while (s_ != s_end_ && !*s_);

      return *this;
    }
    self_type operator++(int junk) { 
      self_type i = *this;
      ++(*this);
      return i;
    }
    reference operator*() const { return *v_; }
    pointer operator->() const { return &(*v_); }
    bool operator==(const self_type& rhs) const { return v_ == rhs.v_; }
    bool operator!=(const self_type& rhs) const { return v_ != rhs.v_; }

  private:
    typename std::vector<T>::const_iterator v_;
    std::vector<bool>::const_iterator s_;
    std::vector<bool>::const_iterator s_end_;
  };

  SubSet(const std::vector<T>& v,
         const std::vector<bool>& s) : v_(v), s_(s) {}
  iterator begin() const {
    return iterator{ v_.begin(), s_.begin(), s_.end() };
  }
  iterator end() const {
    return iterator(v_.end(), s_.end(), s_.end());
  }
private:
  const std::vector<T>& v_;
  const std::vector<bool>& s_;
};
template <typename T>
class PowerSet {
public:
  PowerSet(std::vector<T> v) : v_(std::move(v)) {}

  class iterator {
  public:
    typedef iterator self_type;
    typedef SubSet<T> value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef ptrdiff_t difference_type;

    iterator(const std::vector<T>& v, std::vector<bool> index, bool carry) : v_(v), index_(std::move(index)), carry_(carry) {
    }
    self_type operator++(int) { self_type i = *this; carry_ |= BinaryIncrement(&index_); return i; }
    self_type operator++() {
      carry_ |= BinaryIncrement(&index_);
      return *this; 
    }
    value_type operator*() { return SubSet<T>(v_, index_); }
    bool operator==(const self_type& rhs) { return (carry_ == rhs.carry_) && (index_ == rhs.index_); }
    bool operator!=(const self_type& rhs) { return (carry_ != rhs.carry_) || (index_ != rhs.index_); }
  private:
    const std::vector<T>& v_;
    std::vector<bool> index_;
    bool carry_;
  };
  iterator begin() const {
    return iterator(v_, std::vector<bool>(v_.size(), false), false);
  }
  iterator end() const {
    return iterator(v_, std::vector<bool>(v_.size(), false), true);
  }
private:
  std::vector<T> v_;
};

} // namespace

}

namespace std {
template <>
struct hash<kmap_solver::Group> {
  size_t operator()(const kmap_solver::Group& k) const
  {
    // Compute individual hash values for first,
    // second and third and combine them using XOR
    // and bit shifting:

    return ((hash<unsigned long>()(k.val)
             ^ (hash<unsigned long>()(k.dc) << 1)) >> 1)
      ^ (hash<int>()(k.size) << 1);
  }
};

}

namespace kmap_solver {
std::string Group::SOP(const std::vector<std::string>& variable_names) const
{
  unsigned long d = dc;
  unsigned long v = val;
  std::vector<std::string> terms;
  for (unsigned int i = 0; i < variable_names.size(); i++) {
    // MSB corresponds to largest variable name.
    unsigned int j = variable_names.size() - i - 1;
    if (((1UL << j) & d) == 0) {
      terms.push_back(variable_names[i] + ((
        ((1UL << j) & v) == 0UL) ? "'" : ""));
    }
  }
  return absl::StrJoin(terms, " ");
}

template <typename Iter1, typename Iter2>
class JoinedCollections {
public:
  JoinedCollections(
    Iter1 begin1,
    Iter1 end1,
    Iter2 begin2,
    Iter2 end2) : begin1_(begin1), end1_(end1), begin2_(begin2), end2_(end2) {}
  class iterator {
  public:
    typedef iterator self_type;
    typedef typename Iter1::value_type value_type;
    typedef typename Iter1::reference reference;
    typedef typename Iter1::pointer pointer;
    typedef std::forward_iterator_tag iterator_category;
    typedef ptrdiff_t difference_type;

    iterator(
      Iter1 begin1,
      Iter1 end1,
      Iter2 begin2) : begin1_(begin1), end1_(end1), begin2_(begin2) {}

    self_type operator++(int) { self_type i = *this; ++(*this); return i; }
    self_type operator++() {
      if (begin1_ != end1_) {
        ++begin1_;
      } else {
        ++begin2_;
      }
      return *this;
    }
    reference operator*() const {
      if (begin1_ != end1_) {
        return *begin1_;
      } else {
        return *begin2_;
      }
    }

    pointer operator->() const {
      return &(**this);
    }
    bool operator==(const self_type& rhs) const {
      return begin1_ == rhs.begin1_ &&
        begin2_ == rhs.begin2_;
    }
    bool operator!=(const self_type& rhs) const { return begin1_ != rhs.begin1_ || begin2_ != rhs.begin2_; }

  private:
    Iter1 begin1_;
    Iter1 end1_;
    Iter2 begin2_;
  };

  iterator begin() const {
    return iterator(begin1_, end1_, begin2_);
  }
  iterator end() const {
    return iterator(end1_, end1_, end2_);
  }
private:
  Iter1 begin1_;
  Iter1 end1_;
  Iter2 begin2_;
  Iter2 end2_;
};

template <typename Iterator>
bool GroupCoverage(Iterator grouping_begin, Iterator grouping_end, std::unordered_set<unsigned long> on_vars) {
  // Can be made O(n+m) if I use clever bit hacks.
  return std::all_of(on_vars.begin(), on_vars.end(),
                     [&](unsigned long t)->bool {
    return std::any_of(grouping_begin, grouping_end, [=](const Group& group)->bool {
      return group.Overlap(t);
    }); });
}

template <typename Iterator>
bool RaceHazard(Iterator grouping_begin, Iterator grouping_end) {
  // First get the sets of possible race groups. I make no attempt to de-dup these.
  // There is definetly some way to better filter this so it is not O(n^2), or atleast
  // reduce the cardinality.
  for (Iterator g = grouping_begin; g != grouping_end; ++g) {
    Iterator k = g;
    ++k;
    for (; k != grouping_end; ++k) {
      if (g->RaceHazard(*k) && std::none_of(grouping_begin, grouping_end, [&](
        const Group& d) {
        return d.Overlap(*g) && d.Overlap(*k);
      })) {
        return true;
      }
    }
  }

  // now check if there is no other group that de-races the group.
  return false;
}


template <typename Iterator>
int GroupCost(Iterator grouping_begin, Iterator grouping_end)
{
  int cost = -1;
  for (; grouping_begin != grouping_end; ++grouping_begin) {
    const Group& group = *grouping_begin;
    cost += group.NumGates();
  }
  return cost;
}
template <typename Iterator>
void MinGroup(Iterator starting_grouping_begin,
              Iterator starting_grouping_end,
              std::unordered_set<Group> all_groups, std::unordered_set<unsigned long> on_vars,
              unsigned int* min_cost, std::vector<Group>* min_group) {
  if (!GroupCoverage(starting_grouping_begin,
                     starting_grouping_end,
                     on_vars)) {
    return;
  }
  int cost = GroupCost(starting_grouping_begin, starting_grouping_end);
  if (*min_cost < cost) return;
  if (!RaceHazard(starting_grouping_begin, starting_grouping_end)) {
    min_group->clear();
    for (auto i = starting_grouping_begin; i != starting_grouping_end; ++i) {
      min_group->push_back(*i);
    }
    *min_cost = cost;
    return;
  }
  // if there is a race hazard, find the minimum set without the race (allowing using DC groups
  std::unordered_set<Group> non_group_vals_set = all_groups;
  for (Iterator i = starting_grouping_begin; i != starting_grouping_end; ++i) {
    non_group_vals_set.erase(*i);
  }
  std::vector<Group> non_group_vals(non_group_vals_set.begin(), non_group_vals_set.end());
  for (const auto& set : PowerSet<Group>(non_group_vals)) {
    JoinedCollections<Iterator, decltype(set.begin())> candidate_group(starting_grouping_begin, starting_grouping_end, set.begin(), set.end());
    cost = GroupCost(candidate_group.begin(), candidate_group.end());
    if (cost < *min_cost && !RaceHazard(candidate_group.begin(), candidate_group.end())) {
      min_group->clear();
      for (auto i = starting_grouping_begin; i != starting_grouping_end; ++i) {
        min_group->push_back(*i);
      }
      *min_cost = cost;
    }
  }
  
}

std::unordered_set<Group> CombineAdjacentRectages(
  const std::unordered_set<Group>& last_starts, const std::unordered_set<Group>& last_fins, int num_inputs) {
  std::unordered_set<Group> ret;
  for (auto g = last_starts.begin(); g != last_starts.end(); ++g) {
    for (int j = 0; j < num_inputs; j++) {
      unsigned long bit_mask = 1UL << j;
      // If this bit is already dc, don't need to check.
      if ((g->dc & bit_mask) != 0UL) {
        continue;
      }
      // Make the candidate pair to group with.
      Group candidate_group{ g->val ^ bit_mask, g->dc, num_inputs };
      if (last_fins.find(candidate_group) != last_fins.end()) {
        ret.insert(Group{ g->val & ~bit_mask, g->dc | bit_mask, num_inputs });
      }
    }
  }
  return ret;
}



std::vector<Group> SolveKMap(int num_inputs, std::unordered_set<unsigned long> on_vars, std::unordered_set<unsigned long> dc_vars)
{
  // Only groups with equivalent dont cares, and a hamming distance of 1 can be combined.
  // For size=1, copy over vars with dc = 0.
  std::unordered_set<Group> dc_groups;
  std::unordered_set<Group> c_groups;

  // A group is defined by a pair of ulongs, with 00 being a false, 01 being a true, and
  // 10 being a dc. so 0110, 0001 would correspond to A'BC. 
  for (auto l : on_vars) {
    c_groups.insert(Group{ l, 0, num_inputs });
  }
  for (auto l : dc_vars) {
    dc_groups.insert(Group{ l, 0, num_inputs });
  }
  std::vector<std::unordered_set<Group>> dc_level_groups;
  dc_level_groups.push_back(dc_groups);
  std::vector<std::unordered_set<Group>> c_level_groups;
  c_level_groups.push_back(c_groups);

  unsigned long max_group = (1UL << num_inputs);

  // Iterate by resultant group size.
  for (int i = 1; (1UL << i) <= max_group; i++) {
    dc_level_groups.emplace_back();
    c_level_groups.emplace_back();
    auto  combination = CombineAdjacentRectages(c_level_groups[i - 1], c_level_groups[i - 1], num_inputs);
    c_level_groups.back().insert(combination.begin(), combination.end());
    combination = CombineAdjacentRectages(c_level_groups[i - 1], dc_level_groups[i - 1], num_inputs);
    c_level_groups.back().insert(combination.begin(), combination.end());
    combination = CombineAdjacentRectages(dc_level_groups[i - 1], dc_level_groups[i - 1], num_inputs);
    dc_level_groups.back().insert(combination.begin(), combination.end());;
  }
  // Get all groups in decending order of size.
  std::vector<Group> on_groups;
  for (auto i = c_level_groups.rbegin(); i != c_level_groups.rend(); ++i) {
    on_groups.insert(on_groups.end(), i->begin(), i->end());
  }
  std::unordered_set<Group> deracing_groups;
  // Get all groups that have atleast one elimination.
  for (unsigned int i = 1; i < c_level_groups.size(); i++) {
    deracing_groups.insert(c_level_groups[i].begin(), c_level_groups[i].end());
    deracing_groups.insert(dc_level_groups[i].begin(), dc_level_groups[i].end());
  }


  unsigned int min_cost = std::numeric_limits<unsigned int>::max();
  std::vector<Group> group;
  for (auto grouping : PowerSet<Group>(on_groups)) {
    MinGroup(grouping.begin(),
             grouping.end(),
             deracing_groups, on_vars, &min_cost, &group);
  }
  return group;
}

std::string SOP(const std::vector<Group>& group, const std::vector<std::string>& vars) {
  std::vector<std::string> terms;
  for (const auto& g : group) {
    terms.push_back(g.SOP(vars));
  }
  return absl::StrJoin(terms, " + ");
}
}
