
using BioAlignments
using MetaGraphs
using LightGraphs
using BiDiGraph

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
#infile = "data/example_graph.gfa"
#infile = "/home/cg94xoli/git/MindTheGap_graph/pipeline/genome_graph/failed_simplification/L7Lc20.manually_simplified.gfa"
#outfile = "data/example_out.gfa"



println("Loading graph")
g = readGFA(infile)
writeToGfa(g,outfile*"_read.gfa",kmerSize)  # Should infer overlap

println("Initial graph :")
graph_stats(g)



# Remove all looped gapfillings
g = remove_self_loops!(g)

# Remove reciprocal gapfillings
g = pop_all_bubbles!(g)
writeToGfa(g,outfile*"_bubbles.gfa",kmerSize)  # Should infer overlap

# Merge redundant gapfillings
g = merge_gapfillings!(g,kmerSize)
writeToGfa(g,outfile*"_merged.gfa",kmerSize)  # Should infer overlap

# Detection of linear paths
LinearPaths = findAllLinearPaths(g,kmerSize)
g= merge_all_linear_paths!(g,LinearPaths)

writeToGfa(g,outfile*"_uncut.gfa",kmerSize)  # Should infer overlap


# disconnect branching contigs
v = 1

while v < nv(g)
    global g,v

    if length(outneighbors(g,v)) > 7
        g = cut_edges!(g,v,"R")
    end
    if length(inneighbors(g,v)) > 7
        g = cut_edges!(g,v,"L")
    end
    v = v+1
end


# remove small contigs
v = 1
while v < nv(g)
    # TODO Only if unconnected
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
