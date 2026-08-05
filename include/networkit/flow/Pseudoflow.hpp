/*
 * Pseudoflow.hpp
 *
 *  Created on: 05.08.2026
 *      Author: Roman Erhardt
 * <roman.erhardt@kit.edu>
 */

#ifndef NETWORKIT_FLOW_PSEUDOFLOW_HPP_
#define NETWORKIT_FLOW_PSEUDOFLOW_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>

#include <vector>
#include <iostream>

namespace NetworKit {
    using FlowType = double;

    class Pseudoflow final : public Algorithm
    {

    private:
        class Forest
        {
        public:
            explicit Forest(const index n) : parentLink(n, ParentLink{none, none}),
                                             childLink(n, ChildLink{none, none}),
                                             prevSibling(n, none)
            {
            }

            inline void addEdge(const node parent, const node child, const edgeid edge) noexcept
            {
                parentLink[child].parentEdge = edge;
                parentLink[child].parentVertex = parent;
                prevSibling[child] = none;
                const node fc = childLink[parent].firstChild;
                childLink[child].nextSibling = fc;
                if (fc != none)
                    prevSibling[fc] = child;
                childLink[parent].firstChild = child;
            }

            inline void removeParentEdge(const node vertex) noexcept
            {
                removeChild(vertex);
                parentLink[vertex].parentEdge = none;
                parentLink[vertex].parentVertex = none;
            }

            inline void removeParentOnly(const node vertex) noexcept
            {
                parentLink[vertex].parentEdge = none;
                parentLink[vertex].parentVertex = none;
            }

            inline void changeParentEdge(const node parent, const node child, const edgeid edge) noexcept
            {
                removeChild(child);
                addEdge(parent, child, edge);
            }

            template <typename FUNCTION>
            inline void forAllChildren(const node parent, const FUNCTION &callback) const noexcept
            {
                node child = childLink[parent].firstChild;
                while (child != none)
                {
                    callback(child);
                    child = childLink[child].nextSibling;
                }
            }

            template <typename FUNCTION>
            inline void removeChildren(const node parent, const FUNCTION &remove) noexcept
            {
                node child = childLink[parent].firstChild;
                node prev = none;
                while (child != none)
                {
                    const node next = childLink[child].nextSibling;
                    if (remove(child))
                    {
                        if (prev == none)
                            childLink[parent].firstChild = next;
                        else
                            childLink[prev].nextSibling = next;
                        if (next != none)
                            prevSibling[next] = prev;
                    }
                    else
                    {
                        prev = child;
                    }
                    child = next;
                }
            }

            inline node getParentVertex(const node vertex) const noexcept
            {
                return parentLink[vertex].parentVertex;
            }

            inline node getFirstChild(const node vertex) const noexcept
            {
                return childLink[vertex].firstChild;
            }

            inline edgeid getParentEdge(const node vertex) const noexcept
            {
                return parentLink[vertex].parentEdge;
            }

        private:
            inline void removeChild(const node child) noexcept
            {
                const node ns = childLink[child].nextSibling;
                if (ns != none)
                    prevSibling[ns] = prevSibling[child];
                if (prevSibling[child] == none)
                    childLink[parentLink[child].parentVertex].firstChild = ns;
                else
                    childLink[prevSibling[child]].nextSibling = ns;
            }

            struct ParentLink
            {
                edgeid parentEdge;
                node parentVertex;
            };
            struct ChildLink
            {
                node firstChild;
                node nextSibling;
            };

            std::vector<ParentLink> parentLink;
            std::vector<ChildLink> childLink;
            std::vector<node> prevSibling;
        };

        template <bool CHECK_CONTAINMENT>
        class RootBuckets
        {
            struct Bucket
            {
                std::vector<node> data;
                size_t head = 0;

                bool empty() const noexcept { return head >= data.size(); }

                void push_back(node v) { data.push_back(v); }

                node front() const noexcept { return data[head]; }

                void pop_front() noexcept
                {
                    ++head;

                    // Optional cleanup if the bucket becomes empty
                    if (head == data.size())
                    {
                        data.clear();
                        head = 0;
                    }
                }
            };

        public:
            RootBuckets(const index n) : buckets(n), contained(CHECK_CONTAINMENT ? n : 0, false), maxBucket(-1) {}

            // Should only be called if there are no vertices with distance == buckets.size() - 2.
            // Lift the distances of the components in the highest bucket to n.
            template <typename CALLBACK>
            inline void clearGap(const CALLBACK &callback) noexcept
            {
                while (!buckets[maxBucket].empty())
                {
                    const node root = buckets[maxBucket].front();
                    if constexpr (CHECK_CONTAINMENT)
                        contained[root] = false;
                    buckets[maxBucket].pop_front();
                    callback(root);
                }
                while (maxBucket >= 0 && buckets[maxBucket].empty())
                    maxBucket--;
            }

            // If the only bucket is 0, relabel all the roots to distance 1
            template <typename CALLBACK>
            inline void checkZeroBucket(const CALLBACK &callback) noexcept
            {
                if (maxBucket != 0)
                    return;
                maxBucket = 1;
                while (!buckets[0].empty())
                {
                    const node root = buckets[0].front();
                    buckets[0].pop_front();
                    callback(root);
                    buckets[1].push_back(root);
                }
            }

            inline node extractHighestRoot() noexcept
            {
                if (maxBucket < 0)
                    return none;
                const node root = buckets[maxBucket].front();
                if constexpr (CHECK_CONTAINMENT)
                    contained[root] = false;
                buckets[maxBucket].pop_front();
                while (maxBucket >= 0 && buckets[maxBucket].empty())
                    maxBucket--;
                return root;
            }

            inline void add(const node root, const int dist) noexcept
            {
                if constexpr (CHECK_CONTAINMENT)
                    if (contained[root])
                        return;
                assert(static_cast<index>(dist) < buckets.size());
                maxBucket = std::max(maxBucket, dist);
                buckets[dist].push_back(root);
                if constexpr (CHECK_CONTAINMENT)
                    contained[root] = true;
            }

            inline bool empty() const noexcept
            {
                return maxBucket < 0;
            }

            inline index size() const noexcept
            {
                return static_cast<index>(maxBucket + 1);
            }

        private:
            // Always ensure last bucket is not empty
            // Don't explicitly maintain distance-n bucket
            std::vector<Bucket> buckets;
            std::vector<bool> contained;
            int maxBucket;
        };

    public:
        Pseudoflow(const Graph &G, node s, node t);

        void run() override;

        std::vector<node> getSourceSet() const;

        std::vector<node> getSinkSet() const;

        std::vector<edgeid> getCutEdges() const;

        FlowType getMaxFlow() const;

    private:
        struct Cut
        {
            explicit Cut(const int n, const node source) : inSinkComponent(n, true) { inSinkComponent[source] = false; }

            inline void addToSourceComponent(const node v) noexcept { inSinkComponent[v] = false; }

            std::vector<node> getSourceComponent() const noexcept;

            std::vector<node> getSinkComponent() const noexcept;

            std::vector<bool> inSinkComponent;
        };

        // Saturate source- and sink-incident arcs
        // If this creates excess, make the vertex a strong root
        void initialize() noexcept;

        void runAfterInitialize() noexcept;

        node extractHighestStrongRoot() noexcept;

        // Process component of root in DFS order and set all distances to n
        void liftComponent(const node root) noexcept;

        void processRoot(const node root) noexcept;

        bool tryMerge(const node root, const node u) noexcept;

        std::pair<node, edgeid> findMergerEdge(const node u) noexcept;

        // Invert parent-child relationship along the path from vertex to its root
        void merge(node vertex, node newParentVertex, edgeid newParentEdge) noexcept;

        // Push excess from v towards its root
        void pushExcess(node v, FlowType oldExcess = 1) noexcept;

        void pushFlow(const node from, const node to, const edgeid edge, const edgeid reverseEdge, const FlowType flow) noexcept;

        void incrementDistance(const node vertex) noexcept;

    private:
        const Graph &graph;
        const index n;
        const node sourceVertex;
        const node sinkVertex;
        std::vector<edgeid> reverseEdge;
        std::vector<FlowType> residualCapacity;
        std::vector<int> distance;
        std::vector<int> distanceCount;
        std::vector<FlowType> excess;
        std::vector<index> currentNeighbor;
        Forest forest;
        RootBuckets<true> rootBuckets;
        Cut cut;
        // Reused across processRoot/liftComponent calls to avoid a heap allocation
        // per strong-root DFS (these are called millions of times).
        std::vector<node> dfsStack;
    };
} // namespace NetworKit

#endif // NETWORKIT_FLOW_PSEUDOFLOW_HPP_
