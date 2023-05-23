#include <libsatsuma/Problems/BidirectedGraph.hh>
#include <optional>

namespace Satsuma {

// Which nodes need to be flipped to orient the graph
using Orientation = BidirectedGraph::NodeMap<bool>;

std::unique_ptr<Orientation> try_orient(BidirectedGraph const &bg);

} // namespace Satsuma
