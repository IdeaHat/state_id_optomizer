#include "kmap_solver.h"

#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <unordered_set>
#include <queue>
#include <tuple>
#include <unordered_map>

#include "absl/strings/str_join.h"

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


namespace {


struct AdjacencyListEdge {
  size_t node_id;
  size_t weight;
};



class InOrderPathCusor {
public:
  InOrderPathCusor(std::vector<std::vector<AdjacencyListEdge>> entries, size_t start, size_t end) : edges_(std::move(entries)), start_(start), end_(end), visited_(edges_.size()) {
    priority_queue_.push({ start, 0, {} });
  }

  bool GetNextPath(std::vector<size_t>* nodes, size_t* cost) {
    while (priority_queue_.size()) {
      PriorityQueueEntry entry = priority_queue_.top();
      priority_queue_.pop();
      if (entry.node == end_) {
        *nodes = std::move(entry.path);
        *cost = entry.cost;
        return true;
      }
      entry.path.push_back(entry.node);
      for (const auto& edge : edges_[entry.node]) {
        PriorityQueueEntry new_entry = entry;
        new_entry.node = edge.node_id;
        new_entry.cost += edge.weight;
        priority_queue_.push(new_entry);
      }
    }
    return false;
  }
private:
  struct PriorityQueueEntry {
    size_t node;
    size_t cost;
    std::vector<size_t> path;

    bool operator <(const PriorityQueueEntry& entry) const {
      return entry.cost < cost;
    }
  };
  std::priority_queue<PriorityQueueEntry> priority_queue_;
  std::vector<std::vector<AdjacencyListEdge>> edges_;
  size_t start_;
  size_t end_;
  std::vector<std::unordered_map<size_t, std::unordered_set<size_t>>> visited_;
};

void MakeCostGraphForGroups(
  const std::vector<size_t>& group_counts,
  std::vector<std::vector<AdjacencyListEdge>>* ret_ptr,
  std::vector<size_t>* costs_ptr) {
  // Assume the lowest costs 0, each bin above that costs one more.
  size_t total = 0;
  for (const auto v : group_counts) total += v;
  // What is my cost?
  std::vector<size_t>& costs = *costs_ptr;
  costs = std::vector<size_t>(total, 0);
  for (size_t i = 0, j = 0; i < group_counts.size(); i++) {
    for (size_t k = 0; k < group_counts[i]; k++, j++) {
      costs[j] = i;
    }
  }

  // Let total be the start, total+1 be the finish.
  std::vector<std::vector<AdjacencyListEdge>>& ret = *ret_ptr;
  ret = std::vector<std::vector<AdjacencyListEdge>>(total + 2);

  const size_t start = total;
  const size_t finish = start + 1;
  for (size_t i = 0; i < total; i++) {
    // A 0 connection from start to here.
    ret[start].push_back({ i, 0 });
    // A wight connection from here to finish.
    ret[i].push_back({ finish, costs[i] });
    // A weight+1 connection from here to all previos nodes.
    for (size_t j = 0; j < i; j++) {
      ret[i].push_back({ j, costs[i] + 1 });
    }
  }
}

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
bool GroupingCost(Iterator starting_grouping_begin,
                  Iterator starting_grouping_end,
                  std::unordered_set<unsigned long> on_vars,
                  unsigned int* min_cost, std::vector<Group>* min_group) {
  if (!GroupCoverage(starting_grouping_begin,
                     starting_grouping_end,
                     on_vars)) {
    return false;
  }
  int cost = GroupCost(starting_grouping_begin, starting_grouping_end);

  if (RaceHazard(starting_grouping_begin, starting_grouping_end)) {
    return false;
  }
  min_group->clear();
  for (auto i = starting_grouping_begin; i != starting_grouping_end; ++i) {
    min_group->push_back(*i);
  }
  *min_cost = cost;
  return true;
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


std::tuple<std::vector<std::vector<Group>>, std::vector<std::vector<Group>>> FindAllGroups(
  int num_inputs,
  const std::unordered_set<unsigned long>& on_vars, const std::unordered_set<unsigned long>& dc_vars) {
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


  std::vector<std::unordered_set<Group>> dc_level_groups(num_inputs);
  std::vector<std::unordered_set<Group>> c_level_groups(num_inputs);
  std::vector<size_t> group_counts(num_inputs);

  dc_level_groups.back() = dc_groups;
  c_level_groups.back() = c_groups;
  group_counts.resize(num_inputs, 0);
  group_counts.back() = on_vars.size();

  // Iterate by resultant group size.
  for (int i = num_inputs - 1; i > 0; i--) {
    auto  combination = CombineAdjacentRectages(c_level_groups[i], c_level_groups[i], num_inputs);
    c_level_groups[i - 1].insert(combination.begin(), combination.end());
    combination = CombineAdjacentRectages(c_level_groups[i], dc_level_groups[i], num_inputs);
    c_level_groups[i - 1].insert(combination.begin(), combination.end());
    combination = CombineAdjacentRectages(dc_level_groups[i], dc_level_groups[i], num_inputs);
    dc_level_groups[i - 1].insert(combination.begin(), combination.end());;
    group_counts[i - 1] = c_level_groups[i - 1].size();
  }

  std::tuple<std::vector<std::vector<Group>>, std::vector<std::vector<Group>>> ret;
  for (int i = 0; i < num_inputs; i++) {
    std::get<0>(ret).emplace_back(c_level_groups[i].cbegin(), c_level_groups[i].cend());
    std::get<1>(ret).emplace_back(dc_level_groups[i].cbegin(), dc_level_groups[i].cend());
  }
  return ret;
}

} // namespace

std::string Group::SOP(const std::vector<std::string>& variable_names) const {
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

std::vector<Group> SolveKMap(int num_inputs, const std::unordered_set<unsigned long>& on_vars, const std::unordered_set<unsigned long>& dc_vars)
{

  std::vector<std::vector<Group>> dc_level_groups;
  std::vector<std::vector<Group>> c_level_groups;
  std::tie(c_level_groups, dc_level_groups) = FindAllGroups(num_inputs, on_vars, dc_vars);

  for (int i = 0; i < dc_level_groups.size(); ++i) {
    for (const auto& g : dc_level_groups[i]) {
      c_level_groups[i].push_back(g);
    }

  }

  unsigned int min_cost = std::numeric_limits<unsigned int>::max();
  std::vector<Group> group;

  std::vector<size_t> group_counts;
  for (const auto& v : c_level_groups) {
    group_counts.push_back(v.size());
  }

  std::vector<std::vector<AdjacencyListEdge>> adjacency_list;
  std::vector<size_t> costs;
  MakeCostGraphForGroups(
    group_counts,
    &adjacency_list,
    &costs);
  size_t total = costs.size();
  // flatten the group array.
  std::vector<Group> groups;
  groups.reserve(total);
  for (const auto& l : c_level_groups) {
    for (const auto& g : l) {
      groups.push_back(g);
    }
  }
  auto cursor = InOrderPathCusor(adjacency_list, adjacency_list.size() - 2, adjacency_list.size() - 1);
  std::vector<size_t> nodes;  size_t cost;
  while (cursor.GetNextPath(&nodes, &cost)) {
    std::vector<Group> grouping;
    grouping.reserve(total);
    for (size_t n : nodes) {
      if (n < total) {
        grouping.push_back(groups[n]);
      }
    }

    if (GroupingCost(grouping.begin(),
                     grouping.end(),
                     on_vars, &min_cost, &group)) break;

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
