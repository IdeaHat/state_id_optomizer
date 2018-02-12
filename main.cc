#include <iostream>
#include <vector>
#include <string>
#include <iostream>
#include <string>
#include <fstream>.
#include <cctype>

#include "kmap_solver.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"

int main(int arc, char** argv) {
  using std::string;
  int num_inputs = atoi(argv[1]);
  int num_outs = atoi(argv[2]);
  std::string filename = argv[3];


  std::vector<std::string> input_names;
  std::vector<std::string> output_names;
  std::vector<kmap_solver::StateTransition> transitions;
  int num_states = 0;

  // Read the truth table.

  {
    std::ifstream in(filename);
    std::string line;
    bool first_line = true;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      auto string_values = absl::StrSplit(line, ' ');

      std::vector<string> input_vals;
      std::vector<string> output_vals;
      string from_state;
      string to_state;

      {
        auto i = string_values.begin();
        from_state = string(*i);
        ++i;
        for (int j = 0; j < num_inputs; j++, ++i) {
          input_vals.push_back(string(*i));
        }
        to_state = string(*i);
        ++i;
        for (int j = 0; j < num_outs; j++, ++i) {
          output_vals.push_back(string(*i));
        }
      }

      if (first_line) {
        input_names = std::move(input_vals);
        output_names = std::move(output_vals);
        first_line = false;
        continue;
      } else {
        std::vector<string> ins;
        ins.push_back("0");
        for (int i = 0; i < input_vals.size(); i++) {
          char c = std::tolower(input_vals[i][0]);
          switch (c) {
          case 'x':
          {
            auto i_size = ins.size();
            for (int j = 0; j < i_size; j++) {
              ins.push_back(ins[j]);
              ins[j] += '0';
              ins[j + i_size] += '1';
            }
          }
          break;
          case '0':
          case '1':
            for (int j = 0; j < ins.size(); j++)
              ins[j] += c;
            break;
          default:
            throw "Cannot parse truth table";

          }
        }
        std::vector<size_t> numbers;
        numbers.reserve(ins.size());
        for (const auto& s : ins) {
          numbers.push_back(std::stoi(s, nullptr, 2));
        }
        std::sort(numbers.begin(), numbers.end());
        kmap_solver::Group output;
        output.size = num_outs;
        for (int i = 0; i < num_outs; i++) {
          switch (std::tolower(output_vals[i][0])) {
          case 'x':
            output.dc |= 1 << (num_outs - 1 - i);
            break;
          case '1':
            output.val |= 1 << (num_outs - 1 - i);
            break;
          case '0':
            break;
          default:
            throw "Cannot parse truth table";
          }
        }
        int  from_state_num = std::stoi(from_state);
        int  to_state_num =
          std::tolower(to_state[0]) == 'x' ?
          -1 : std::stoi(to_state);
        num_states = std::max(num_states,
                              std::max(from_state_num, to_state_num));
        for (auto n : numbers) {
          transitions.push_back(
            { from_state_num,
            n,
            to_state_num,
            output });
        }
      }
    }
  }
  using namespace kmap_solver;
  std::vector<size_t> state_ids;
  std::vector<std::vector<std::vector<Group>>> flip_flop_inputs;
  std::vector<std::vector<Group>> outputs;
  const auto flip_flip = SNotR();
  // Solve the KMaps.
  AssignStates(
    flip_flip,
    num_states + 1,
    num_inputs,
    num_outs,
    transitions,
    &state_ids,
    &flip_flop_inputs,
    &outputs);

  PrintStateAssignmentSolution(
    state_ids,
    flip_flop_inputs,
    outputs,
    input_names,
    output_names,
    flip_flip);

  return 0;
}