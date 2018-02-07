#include <iostream>
#include <vector>
#include <string>

#include "kmap_solver.h"
#include "abseil-cpp/absl/strings/str_join.h"

int main(int arc, char** argv) {
  std::cout << absl::StrJoin(std::vector<std::string>{kmap_solver::Foo(), "absl_test"}, "+") << std::endl;
  return 0;
}