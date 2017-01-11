#include <gtest/gtest.h>
#include "wmrde/algebra/matrix.h"

//using namespace wmrde;

// Declare a test
TEST(TestSuite, addTwoInts)
{
  EXPECT_EQ(wmrde::addTwoInts(2,2),4);
}

// Declare another test
TEST(TestSuite, addTwoIntsFail)
{
  EXPECT_EQ(wmrde::addTwoInts(2,2),5);
}

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
