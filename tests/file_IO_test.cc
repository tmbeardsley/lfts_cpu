#include "../src/file_IO.h"
#include <gtest/gtest.h>


// Tests for file_IO::isValidFile()
TEST(fileIOTests, isValidFile) {

  // Should return false for empty string
  EXPECT_FALSE(file_IO::isValidFile(""));

  // Should return false for a fake directory
  EXPECT_FALSE(file_IO::isValidFile("./a/fake/dir/nonExistent.txt"));

  // Should return true for an existing file
  EXPECT_TRUE(file_IO::isValidFile("../src/file_IO.h"));
}