#include <libsatsuma/Solvers/OrientBinet.hh>
#include <lemon/dfs.h>
#include <memory>

namespace Satsuma {

namespace detail {

class not_orientable_error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

template<typename GraphT>
class OrientVisitor : public lemon::DfsVisitor<GraphT>
{
public:
    using Arc = typename GraphT::Arc;
    using Edge = typename GraphT::Edge;
    OrientVisitor(BidirectedGraph const& _bg,
            Orientation &_ori)
        : bg_(_bg)
        , ori_{_ori}
    {}
    void discover(Arc const &a)
    {
        auto e = Edge{a};
        auto src = bg_.g.source(a);
        auto tgt = bg_.g.target(a);
#if 0
        if (src == bg_.outer_node
         || tgt == bg_.outer_node) {
            return;
        }
#endif

        bool is_bidir = (bg_.u_head[a] == bg_.v_head[a]);
        ori_[tgt] = ori_[src] ^ is_bidir;
#if 0
        std::cout << "disco ";
        debug_print(a);
#endif
    }
    void examine(Arc const &a)
    {
        auto e = Edge{a};
        auto src = bg_.g.source(a);
        auto tgt = bg_.g.target(a);
#if 0
        if (src == bg_.outer_node
         || tgt == bg_.outer_node) {
            return;
        }
#endif
        bool is_bidir = (bg_.u_head[a] == bg_.v_head[a]);
        if ((ori_[tgt] ^ ori_[src]) != is_bidir) {
            // TODO: there should be a better way to abort.
            std::cout << " exam ";
            debug_print(a);
            throw not_orientable_error("");
        };
    }
private:
    void debug_print(Arc const &a)
    {
        auto e = Edge{a};
        auto src = bg_.g.source(a);
        auto tgt = bg_.g.target(a);
        bool is_fw = (bg_.g.source(a) == bg_.g.u(e));
        const auto &src_head = is_fw ? bg_.u_head : bg_.v_head;
        const auto &tgt_head = is_fw ? bg_.v_head : bg_.u_head;
        std::cout << "arc " << bg_.g.id(a)
            <<  ": " << bg_.g.id(bg_.g.source(a)) << " (" << src_head[a] << ")"
            << " - " << bg_.g.id(bg_.g.target(a)) << " (" << tgt_head[a] << ")"
            << ", ori " << ori_[src]
            << " - " << ori_[tgt]
            << std::endl;
    }
    BidirectedGraph const& bg_;
    Orientation &ori_;
};
} // namespace detail


std::unique_ptr<Orientation> try_orient(BidirectedGraph const &bg)
{
    auto ori = std::make_unique<Orientation>(bg.g, true);
    auto &g = bg.g;
    auto vis = detail::OrientVisitor<BidirectedGraph::GraphT>{bg, *ori};
    auto dfs = lemon::DfsVisit<BidirectedGraph::GraphT, decltype(vis)>(bg.g, vis);
    try {
        dfs.run();
    } catch (detail::not_orientable_error const&) {
        return {};
    }
    return ori;
}

} // namespace Satsuma
