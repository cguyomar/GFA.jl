
using BioAlignments
using MetaGraphs

include("Utils.jl")
include("GraphUtils.jl")
include("gfa_IO.jl")
include("Bubbles.jl")
include("Path.jl")

g = readGFA("data/example_graph.gfa")

v=1
while v < nv(g)
# while v < 34
    nodeName = get_prop(g,v,:name)
    print(v)
    print(" ")
    g = pop_bubble!(g,v)
    if nodeName==get_prop(g,v,:name)
        v = v+1
    end
end
writeToGfa(g,"data/nodup.gfa",63) # Should infer overlap

# Detection of linear paths
LinearPaths = findAllLinearPaths(g)
for path in LinearPaths
    merge_path!(g,path)
end

writeToGfa(g,"data/output.gfa",63) # Should infer overlap
