#ifndef KMAP_SOLVER_H
#define KMAP_SOLVER_H

#include <string>
#include <unordered_set>

namespace kmap_solver {

struct Group {
  unsigned long val;
  unsigned long dc;
  int size;
  unsigned long bit_mask() const {
    return (1UL << size) - 1UL;
  }

  unsigned long OverlapLong(const Group& t) const {
    return ((t.val | val) & (~dc) & (~t.dc) & ((~t.val) | (~val))) & bit_mask();
  }
  bool Overlap(const Group& t) const {
    return OverlapLong(t) == 0;
  }
  bool Overlap(unsigned long v) const
  {
    return Overlap(Group{ v, 0UL, size });
  }
  bool RaceHazard(const Group& t) const {
    // Are there any races between two groups? You can tell there is a race if the two groups are disjoint
    // by only one term.
    unsigned long l = OverlapLong(t);
    return l != 0 && ((l & (l - 1)) == 0);
  }
  int NumDC() const {
    int count = 0;
    for (int i = 0; i < size; i++) {
      if ((dc & (1UL << i)) != 0) {
        count++;
      }
    }
    return count;
  }

  int NumGates() const {
    return size - NumDC() - 1;
  }

  bool operator ==(const Group& other) const {
    return other.val == val && other.dc == dc && other.size == size;
  }
  std::string SOP(const std::vector<std::string>& variable_names) const;
};

std::vector<Group> SolveKMap(int num_inputs, const std::unordered_set<unsigned long>& on_vars, const std::unordered_set<unsigned long>& dc_vars);

std::string SOP(const std::vector<Group>& group, const std::vector<std::string>& vars);
}
#endif