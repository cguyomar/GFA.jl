
using BioAlignments
using MetaGraphs

include("Utils.jl")
include("GraphUtils.jl")
include("gfa_IO.jl")
include("Bubbles.jl")
include("Path.jl")
include("redundant_gapfillings.jl")

println("Loading graph")
g = readGFA("data/minia_k71_abundancemin_2_filtered_300_gapfilling_k61_abundancemin_auto.gfa")
kmerSize = 63

println("Initial graph :")
graph_stats(g)

# Remove all looped gapfillings
v=1
while v < nv(g)
    if length(intersect(collect(keys(neighbors(g,v,"R"))),collect(keys(neighbors(g,v,"L"))))) >0
        rem_vertex!(g,v)
    else v=v+1
    end
end

# Remove reciprocal gapfillings
v=1
while v < nv(g)
    nodeName = get_prop(g,v,:name)
    print(v)
    print(" ")
    g = pop_bubble!(g,v)
    if nodeName==get_prop(g,v,:name)
        v = v+1
    end
end

# Merge redundant gapfillings
v = 1
while v < nv(g)
    merge_redundant_gapfillings!(g,v,"L")
    merge_redundant_gapfillings!(g,v,"R")
    v = v+1
end

# Detection of linear paths
LinearPaths = findAllLinearPaths(g)
i=1
for path in LinearPaths
    print(i)
    merge_path!(g,path)
    i=i+1
end

# disconnect branching contigs
v = 1
while v < nv(g)
    if length(neighbors(g,v,"R")) >= 5
        g = cut_edges!(g,v,"R")
    end
    if length(neighbors(g,v,"R")) >= 5
        g = cut_edges!(g,v,"R")
    end
    v = v+1
end

# remove small contigs
v = 1
while v < nv(g)
    if length(get_prop(g,v,:seq)) < 100
        rem_vertex!(g,v)
    else
        v = v+1
    end
end

println("Cleaned graph :")
graph_stats(g)

writeToGfa(g,"data/output.gfa",63)  # Should infer overlap
