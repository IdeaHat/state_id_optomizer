#ifndef KMAP_SOLVER_H
#define KMAP_SOLVER_H

#include <string>
#include <unordered_set>

namespace kmap_solver {

struct Group {
  size_t val = 0;
  size_t dc = 0;
  int size = 0;
  size_t bit_mask() const {
    return (1UL << size) - 1UL;
  }

  size_t OverlapLong(const Group& t) const {
    return ((t.val | val) & (~dc) & (~t.dc) & ((~t.val) | (~val))) & bit_mask();
  }
  bool Overlap(const Group& t) const {
    return OverlapLong(t) == 0;
  }
  bool Overlap(size_t v) const
  {
    return Overlap(Group{ v, 0UL, size });
  }
  bool RaceHazard(const Group& t) const {
    // Are there any races between two groups? You can tell there is a race if the two groups are disjoint
    // by only one term.
    size_t l = OverlapLong(t);
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

std::vector<Group> SolveKMap(int num_inputs, const std::unordered_set<size_t>& on_vars, const std::unordered_set<size_t>& dc_vars,
                             size_t* cost_ptr = nullptr);

std::string SOP(const std::vector<Group>& group, const std::vector<std::string>& vars);

struct StateTransition {
  // Must be above 0.
  ptrdiff_t from_state;
  // Inputs.
  size_t inputs;
  // To state. Less than 0 indicates don't care.
  ptrdiff_t to_state;
  // Outputs.
  Group outputs;
};

struct FlipFlop {
  size_t gate_count;
  size_t input_count;
  Group transitions[4];
};

inline FlipFlop SNotR() {
  return { 2, 2,
  {{0b00, 0b01, 2},
  {0b11, 0b00, 2},
  {0b00, 0b00, 2},
  {0b01, 0b10, 2}} };
}

inline size_t transition_index(const bool from_state, const bool to_state) {
  return (from_state ? 0b10 : 0b00) | (to_state ? 0b01 : 0b00);
}

void AssignStates(
  const FlipFlop& flip_flop,
  size_t num_states,
  size_t num_inputs,
  size_t num_outputs,
  const std::vector<StateTransition>& state_transitions,
  std::vector<size_t>* state_ids,
  std::vector<std::vector<std::vector<Group>>>* flip_flop_inputs,
  std::vector<std::vector<Group>>* outputs);

}
#endif