
using BioAlignments
using MetaGraphs

include("Utils.jl")
include("GraphUtils.jl")
include("gfa_IO.jl")
include("Bubbles.jl")
include("Path.jl")
include("redundant_gapfillings.jl")
include("graph2contig.jl")

if length(ARGS)!=3
    println("Usage : ")
    println("julia GFA_analysis.jl infile outfile kmerSize")
end

infile = ARGS[1]
outfile = ARGS[2]
kmerSize = parse(Int,ARGS[3])


println("Loading graph")
g = readGFA(infile)

println("Initial graph :")
graph_stats(g)



# Remove all looped gapfillings
g = remove_self_loops!(g)
# Remove reciprocal gapfillings
g = pop_all_bubbles!(g)

# Merge redundant gapfillings
g = merge_gapfillings!(g,kmerSize)

# Detection of linear paths
LinearPaths = findAllLinearPaths(g,kmerSize)
g= merge_all_linear_paths!(g,LinearPaths)

writeToGfa(g,outfile*"_uncut.gfa",kmerSize)  # Should infer overlap


# disconnect branching contigs
v = 1

while v < nv(g)
    global g,v

    if length(neighbors(g,v,"R")) > 7
        g = cut_edges!(g,v,"R")
    end
    if length(neighbors(g,v,"L")) > 7
        g = cut_edges!(g,v,"L")
    end
    v = v+1
end


# remove small contigs
v = 1
while v < nv(g)
    global g,v
    if length(get_prop(g,v,:seq)) < 100
        rem_vertex!(g,v)
    else
        v = v+1
    end
end

println("Cleaned graph :")
graph_stats(g)

writeToGfa(g,outfile,kmerSize)  # Should infer overlap

graph2contig(g,dirname(outfile),kmerSize)
