
using BioAlignments
using MetaGraphs

include("Utils.jl")
include("GraphUtils.jl")
include("gfa_IO.jl")
include("Bubbles.jl")

g = readGFA("data/example_graph.gfa")

v=1
while v < nv(g)
    print(v)
    print(" ")
    g = pop_bubble!(g,v)
    v = v+1
end
writeToGfa(g,"data/output.gfa",63) # Should infer overlap
