/*
 * Pseudoflow.hpp
 *
 *  Created on: 05.08.2026
 *      Author: Roman Erhardt
 * <roman.erhardt@kit.edu>
 */

#include <stdexcept>

#include <networkit/flow/Pseudoflow.hpp>

namespace NetworKit {

NetworKit::Pseudoflow::Pseudoflow(const Graph &G, node s, node t)
    : graph(G), n(G.upperNodeIdBound()), sourceVertex(s), sinkVertex(t), distance(n, 0),
      distanceCount(n, 0), excess(n, 0), currentNeighbor(n, 0), forest(n), rootBuckets(n),
      cut(n, sourceVertex) {
    if (!graph.hasNode(s) || !graph.hasNode(t))
        throw std::runtime_error("Pseudoflow requires source and sink nodes in the graph!");
    if (s == t)
        throw std::runtime_error("Pseudoflow requires source and sink nodes to be different!");
    if (!graph.isDirected())
        throw std::runtime_error("Pseudoflow requires a directed graph!");
    if (!graph.hasEdgeIds())
        throw std::runtime_error("Pseudoflow requires indexed edges!");

    const index m = graph.upperEdgeIdBound();
    reverseEdge.resize(m, none);
    residualCapacity.resize(m);
    graph.forEdges([&](node u, node v) {
        if (graph.weight(u, v) < 0.0)
            throw std::runtime_error("Pseudoflow requires non-negative capacities!");
        edgeid edge = graph.edgeId(u, v);
        if (!graph.hasEdge(v, u))
            throw std::runtime_error("Pseudoflow requires reverse residual arcs!");
        reverseEdge[edge] = graph.edgeId(v, u);
        residualCapacity[edge] = graph.isWeighted() ? graph.weight(u, v) : 1.0;
    });

    dfsStack.reserve(n);
}

void NetworKit::Pseudoflow::run() {
    hasRun = false;
    std::fill(distance.begin(), distance.end(), 0);
    std::fill(distanceCount.begin(), distanceCount.end(), 0);
    std::fill(excess.begin(), excess.end(), 0.0);
    std::fill(currentNeighbor.begin(), currentNeighbor.end(), 0);
    forest = Forest(n);
    rootBuckets = RootBuckets<true>(n);
    cut = Cut(n, sourceVertex);
    graph.forEdges([&](node u, node v) {
        residualCapacity[graph.edgeId(u, v)] = graph.isWeighted() ? graph.weight(u, v) : 1.0;
    });
    initialize();
    runAfterInitialize();
    hasRun = true;
}

std::vector<edgeid> Pseudoflow::getCutEdges() const {
    assureFinished();
    std::vector<edgeid> edges;
    graph.forNodes([&](node u) {
        if (cut.inSinkComponent[u])
            return;

        graph.forEdgesOf(u, [&](node u, node v, edgeid edge) {
            if (!cut.inSinkComponent[v])
                return;

            // assert(residualCapacity[edge] == 0);
            edges.emplace_back(edge);
        });
    });
    return edges;
}

std::vector<node> Pseudoflow::getSourceSet() const {
    assureFinished();
    std::vector<node> result;
    graph.forNodes([&](node u) {
        if (!cut.inSinkComponent[u])
            result.push_back(u);
    });
    return result;
}

std::vector<node> Pseudoflow::getSinkSet() const {
    assureFinished();
    std::vector<node> result;
    graph.forNodes([&](node u) {
        if (cut.inSinkComponent[u])
            result.push_back(u);
    });
    return result;
}

FlowType NetworKit::Pseudoflow::getMaxFlow() const {
    assureFinished();
    FlowType flow = 0;

    for (const node u : graph.nodeRange()) {
        if (cut.inSinkComponent[u])
            continue;

        graph.forEdgesOf(u, [&](node v) {
            if (!cut.inSinkComponent[v])
                return;
            flow += graph.weight(u, v);
        });
    }
    return flow;
}

std::vector<node> Pseudoflow::Cut::getSourceComponent() const noexcept

{
    std::vector<node> component;
    for (index i = 0; i < inSinkComponent.size(); i++) {
        if (!inSinkComponent[i])
            component.emplace_back(i);
    }
    return component;
}

std::vector<node> Pseudoflow::Cut::getSinkComponent() const noexcept

{
    std::vector<node> component;
    for (index i = 0; i < inSinkComponent.size(); i++) {
        if (inSinkComponent[i])
            component.emplace_back(i);
    }
    return component;
}

void Pseudoflow::initialize() noexcept {
    graph.forEdgesOf(sourceVertex, [&](node u, node v, edgeid edge) {
        const edgeid reverse_edge = reverseEdge[edge];
        const FlowType capacity = graph.weight(sourceVertex, v);
        excess[v] += capacity;
        residualCapacity[edge] = 0;
        residualCapacity[reverse_edge] += capacity;
    });
    graph.forEdgesOf(sinkVertex, [&](node u, node v, edgeid edge) {
        const edgeid reverse_edge = reverseEdge[edge];
        const FlowType capacity = graph.weight(v, sinkVertex);
        excess[v] -= capacity;
        residualCapacity[reverse_edge] = 0;
        residualCapacity[edge] += capacity;
    });

    for (node u = 0; u < n; ++u) {
        if (u == sourceVertex || u == sinkVertex)
            continue;
        if (excess[u] >= 0) {
            distance[u] = 1;
            distanceCount[1]++;
            rootBuckets.add(u, 1);
        }
    }
    distanceCount[0] = n - 2 - distanceCount[1];
}

void Pseudoflow::runAfterInitialize() noexcept {
    node strongRoot = extractHighestStrongRoot();
    while (strongRoot != none) {
        processRoot(strongRoot);
        strongRoot = extractHighestStrongRoot();
    }
}

node Pseudoflow::extractHighestStrongRoot() noexcept {
    while (rootBuckets.size() > 1 && distanceCount[rootBuckets.size() - 2] == 0) {
        rootBuckets.clearGap([&](const node root) { liftComponent(root); });
    }
    rootBuckets.checkZeroBucket([&](const node root) { incrementDistance(root); });
    return rootBuckets.extractHighestRoot();
}

void Pseudoflow::liftComponent(const node root) noexcept {
    auto &stack = dfsStack;
    stack.clear();
    stack.push_back(root);
    while (!stack.empty()) {
        const node v = stack.back();
        stack.pop_back();
        distanceCount[distance[v]]--;
        distance[v] = n;
        cut.addToSourceComponent(v);
        forest.forAllChildren(v, [&](const node child) { stack.push_back(child); });
    }
}

void Pseudoflow::processRoot(const node root) noexcept {
    // Do DFS among subtree that has the same distance as the root
    // If a vertex has an outgoing admissible edge, merge
    // Else, increment the distance
    const int rootDistance = distance[root];
    auto &stack = dfsStack;
    stack.clear();
    stack.push_back(root);
    while (!stack.empty()) {
        const node v = stack.back();
        stack.pop_back();
        if (tryMerge(root, v))
            return;
        bool hasChildren = false;
        forest.forAllChildren(v, [&](const node child) {
            if (distance[child] == rootDistance) {
                stack.push_back(child);
                hasChildren = true;
            }
        });
        if (!hasChildren) {
            incrementDistance(v);
        }
    }
    rootBuckets.add(root, distance[root]);
}

bool Pseudoflow::tryMerge(const node root, const node u) noexcept {
    auto [v, e] = findMergerEdge(u);
    if (v == none)
        return false;
    merge(u, v, e);
    pushExcess(root);
    return true;
}

std::pair<node, edgeid> Pseudoflow::findMergerEdge(const node u) noexcept {
    const int targetDist = distance[u] - 1;

    for (std::pair<node, edgeid> neigh = graph.getIthNeighborWithId(u, currentNeighbor[u]);
         currentNeighbor[u] < graph.degree(u);
         neigh = graph.getIthNeighborWithId(u, ++currentNeighbor[u])) {
        const node &v = neigh.first;
        const edgeid &edge = neigh.second;
        if (residualCapacity[edge] == 0)
            continue;
        if (distance[v] != targetDist)
            continue;
        if (v == sourceVertex || v == sinkVertex)
            continue;
        return std::make_pair(v, edge);
    }
    return std::make_pair(none, none);
}

void Pseudoflow::merge(node vertex, node newParentVertex, edgeid newParentEdge) noexcept {
    while (forest.getParentVertex(vertex) != none) {
        const node oldParentVertex = forest.getParentVertex(vertex);
        const edgeid oldParentEdge = reverseEdge[forest.getParentEdge(vertex)];
        forest.changeParentEdge(newParentVertex, vertex, newParentEdge);
        newParentVertex = vertex;
        vertex = oldParentVertex;
        newParentEdge = oldParentEdge;
    }
    forest.addEdge(newParentVertex, vertex, newParentEdge);
}

void Pseudoflow::pushExcess(node v, FlowType oldExcess) noexcept {
    while (excess[v] > 0 && forest.getParentVertex(v) != none) {
        const edgeid edge = forest.getParentEdge(v);
        const edgeid reverse_edge = reverseEdge[edge];
        const node parent = forest.getParentVertex(v);
        oldExcess = excess[parent];
        if (residualCapacity[edge] >= excess[v]) {
            pushFlow(v, parent, edge, reverse_edge, excess[v]);
        } else {
            pushFlow(v, parent, edge, reverse_edge, residualCapacity[edge]);
            // Parent edge is saturated, so turn vertex into a strong root
            forest.removeParentEdge(v);
            rootBuckets.add(v, distance[v]);
        }
        v = parent;
    }
    if (excess[v] >= 0 && oldExcess < 0) {
        // Weak root has become strong
        rootBuckets.add(v, distance[v]);
    }
}

void Pseudoflow::pushFlow(const node from, const node to, const edgeid edge,
                          const edgeid reverseEdge, const FlowType flow) noexcept {
    residualCapacity[edge] -= flow;
    residualCapacity[reverseEdge] += flow;
    excess[from] -= flow;
    excess[to] += flow;
}

void Pseudoflow::incrementDistance(const node vertex) noexcept {
    distanceCount[distance[vertex]]--;
    distance[vertex]++;
    distanceCount[distance[vertex]]++;
    currentNeighbor[vertex] = 0;
}
} // namespace NetworKit
