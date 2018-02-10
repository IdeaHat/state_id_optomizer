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
     "t' v'" ,
     "t v" ,
     "d" )) << sop;
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

TEST(KMapSolveTest, FourVarsRaceWithDC) {

  auto groups = SolveKMap(4,
                          std::unordered_set<unsigned long>{ 0b0000, 0b0001, 0b1111, 0b1110 },
                          std::unordered_set<unsigned long>{ 0b0100, 0b0101, 0b0111, 0b0110});

  std::string sop = SOP(groups, { "a", "b", "c", "d" });

  EXPECT_THAT(
    groups,
    testing::UnorderedElementsAre(
      Group{ 0b0000, 0b0101, 4 },
      Group{ 0b0110, 0b0110, 4 },
      Group{ 0b0100, 0b0011, 4 })) << sop;
  // split the terms by + and trim.
  std::vector<std::string> terms = absl::StrSplit(sop, "+");
  for (std::string& t : terms) {
    t = trim(t);
  }

  EXPECT_THAT(terms, testing::UnorderedElementsAre(
    "a' c'",
    "b c",
    "a' b")) << sop;
}

TEST(KMapSolveTest, FiveVarsRaceWithDC) {

  auto groups = SolveKMap(5,
                          std::unordered_set<unsigned long>{ 0b00000, 0b00010, 0b01000, 0b01111, 0b11011, 0b11111 },
                          std::unordered_set<unsigned long>{ 0b01010, 0b01011});

  std::string sop = SOP(groups, { "a", "b", "c", "d", "e" });

  EXPECT_THAT(
    groups,
    testing::UnorderedElementsAre(
      Group{ 0b00000, 0b01010, 5 },
      Group{ 0b01011, 0b10100, 5 },
      Group{ 0b01010, 0b00001, 5 })) << sop;
  // split the terms by + and trim.
  std::vector<std::string> terms = absl::StrSplit(sop, "+");
  for (std::string& t : terms) {
    t = trim(t);
  }

  EXPECT_THAT(terms, testing::UnorderedElementsAre(
    "a' c' e'",
    "b d e",
    "a' b c' d")) << sop;
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
    "b' c e' f")) << sop;
}


}
}