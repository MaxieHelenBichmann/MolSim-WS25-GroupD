#include "io/checkpointWriter/YAMLWriterCP.h"

#include <gtest/gtest.h>

namespace mol_sim {
class YAMLWriterCPTest : public testing::Test {
   protected:
    YAMLWriterCP writer;
};

// first some functions testing that the writer writes valid YAML files for tthat project, with valid syntax for
// settings and particles
TEST_F(YAMLWriterCPTest, testWrite) {}

// then some tests that read back the written checkpoint files and check for consistency of data of the whole
// simulation state
TEST_F(YAMLWriterCPTest, testConsistentState) {}
}  // namespace mol_sim