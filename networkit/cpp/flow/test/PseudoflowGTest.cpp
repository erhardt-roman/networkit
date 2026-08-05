/* PseudoflowGTest.cpp */

#include <gtest/gtest.h>

#include <networkit/flow/Pseudoflow.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

TEST(PseudoflowGTest, ComputesUnitFlowWithExplicitResidualArc) {
    Graph graph(2, true, true);
    graph.addEdge(0, 1, 1.0);
    graph.addEdge(1, 0, 0.0);
    graph.indexEdges();

    Pseudoflow algorithm(graph, 0, 1);
    EXPECT_THROW(algorithm.getMaxFlow(), std::runtime_error);

    algorithm.run();

    EXPECT_DOUBLE_EQ(algorithm.getMaxFlow(), 1.0);
    EXPECT_EQ(algorithm.getSourceSet(), (std::vector<node>{0}));
    EXPECT_EQ(algorithm.getSinkSet(), (std::vector<node>{1}));
    EXPECT_EQ(algorithm.getCutEdges(), (std::vector<edgeid>{graph.edgeId(0, 1)}));

    algorithm.run();
    EXPECT_DOUBLE_EQ(algorithm.getMaxFlow(), 1.0);
}

} // namespace NetworKit
