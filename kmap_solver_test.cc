#include "kmap_solver.h"

#include <unordered_set>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "absl/strings/str_split.h"

namespace kmap_solver {
namespace {
std::string trim(std::string& str) {
  size_t first = str.find_first_not_of(' ');
  size_t last = str.find_last_not_of(' ');
  return str.substr(first, (last - first + 1));
}

TEST(KMapSolveTest, ThreeVarsNoRace) {
  /*
    t v d q
    0 0 0 1
    0 0 1 1
    0 1 0 0
    0 1 1 x
    1 0 0 0
    1 0 1 1
    1 1 0 1
    1 1 1 x
    */
  auto groups = SolveKMap(3, std::unordered_set<unsigned long>{ 0b000, 0b001, 0b101, 0b110 },
                          std::unordered_set<unsigned long>{ 0b011, 0b111 });
  EXPECT_THAT(
    groups,
    testing::UnorderedElementsAre(
      Group{ 0b000, 0b001, 3 },
      Group{ 0b110, 0b001, 3 },
      Group{ 0b001, 0b110, 3 }));

  std::string sop = SOP(groups, { "t", "v", "d" });

  // split the terms by + and trim.
  std::vector<std::string> terms = absl::StrSplit(sop, "+");
  for (std::string& t : terms) {
    t = trim(t);
  }

  EXPECT_THAT(terms, testing::UnorderedElementsAre(
    "t' v'",
    "t v",
    "d")) << sop;
}

TEST(KMapSolveTest, ThreeVarsRace) {
  /*
  a b c q
  0 0 0 1
  0 0 1 1
  0 1 0 0
  0 1 1 1
  1 0 0 0
  1 0 1 0
  1 1 0 0
  1 1 1 1
  */
  auto groups = SolveKMap(3, std::unordered_set<unsigned long>{ 0b000, 0b001, 0b111, 0b011 }, {});

  std::string sop = SOP(groups, { "a", "b", "c" });

  EXPECT_THAT(
    groups,
    testing::UnorderedElementsAre(
      Group{ 0b000, 0b001, 3 },
      Group{ 0b011, 0b100, 3 },
      Group{ 0b001, 0b010, 3 })) << sop;
  // split the terms by + and trim.
  std::vector<std::string> terms = absl::StrSplit(sop, "+");
  for (std::string& t : terms) {
    t = trim(t);
  }

  EXPECT_THAT(terms, testing::UnorderedElementsAre(
    "a' b'",
    "b c",
    "a' c")) << sop;
}


TEST(KMapSolveTest, SixVarsRaceWithDC) {

  auto groups = SolveKMap(6,
                          std::unordered_set<unsigned long>{
    0b000000, 0b000001, 0b001000,
      0b001111, 0b011101, 0b011111,
      0b100001, 0b100011, 0b101011,
      0b101100, 0b111100, 0b111101},

                          std::unordered_set<unsigned long>{ 0b001001, 0b001101, 0b101001, 0b101101 });

  std::string sop = SOP(groups, { "a", "b", "c", "d", "e", "f" });

  EXPECT_THAT(
    groups,
    testing::UnorderedElementsAre(
      Group{ 0b000000, 0b001001, 6 },
      Group{ 0b001101, 0b010010, 6 },
      Group{ 0b100001, 0b001010, 6 },
      Group{ 0b101100, 0b010001, 6 },
      Group{ 0b001001, 0b100100, 6 }
  )) << sop;
  // split the terms by + and trim.
  std::vector<std::string> terms = absl::StrSplit(sop, "+");
  for (std::string& t : terms) {
    t = trim(t);
  }

  EXPECT_THAT(terms, testing::UnorderedElementsAre(
    "a' b' d' e'",
    "a' c d f",
    "a c d e'",
    "b' c e' f",
    "a b' d' f")) << sop;
}

TEST(KMapSolveTest, ReallyExpensive6var) {

  std::unordered_set<unsigned long> ones{
    0b000000, 0b000011, 0b001001, 0b001010, 0b011000, 0b011010, 0b010001, 0b010010,
      0b000101, 0b000110, 0b001100, 0b001111, 0b011101, 0b011110, 0b010100, 0b010111,
      0b100100, 0b100111, 0b101101, 0b101110, 0b111100, 0b111110, 0b110101, 0b110110,
      0b100001, 0b100010, 0b101000, 0b101011, 0b111001, 0b111010, 0b110000, 0b110011 };
  auto terms = SolveKMap(6,
                         ones,
                         std::unordered_set<unsigned long>{ });

  std::vector<Group> groups;
  for (auto l : ones) {
    groups.push_back({ l, 0, 6 });
  }
  EXPECT_THAT(terms, testing::UnorderedElementsAreArray(groups));
}


}
}