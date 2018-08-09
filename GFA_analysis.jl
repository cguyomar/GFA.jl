
using BioAlignments
using MetaGraphs
using BioSequences

include("Utils.jl")
include("GraphUtils.jl")
include("gfa_IO.jl")
include("Bubbles.jl")
include("Path.jl")
include("redundant_gapfillings.jl")

if length(ARGS)!=3
    println("Usage : ")
    println("julia GFA_analysis.jl infile outfile kmerSize")
end

infile = ARGS[1]
outfile = ARGS[2]
kmerSize = parse(Int,ARGS[3])


println("Loading graph")
infile = "data/ArPo28.gfa"
g = readGFA(infile)

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
LinearPaths = findAllLinearPaths(g,kmerSize)
i=1
for path in LinearPaths
    merge_path!(g,path)
    i=i+1
end

# disconnect branching contigs
v = 1
while v < nv(g)
    if length(neighbors(g,v,"R")) >= 5
        g = cut_edges!(g,v,"R")
    end
    if length(neighbors(g,v,"L")) >= 5
        g = cut_edges!(g,v,"L")
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

writeToGfa(g,outfile,kmerSize)  # Should infer overlap



# output best sequences
connectedComponents = weakly_connected_components(g.graph)
println("Number of connected components : $(length(connectedComponents))")
nbComponent = 0
for component in connectedComponents
    nbComponent += 1
    found = false
    paths = Vector{Path}()
    for node in component
        if isDeadEnd(g,node)==false
            continue
        else
            found = true
            paths = vcat(paths,find_all_paths(g,node,isDeadEnd(g,node)))
        end
    end
    if found == false # No deadend found -> circle
        paths=[find_all_paths(g,component[1],"+")]
    end
    paths = remove_duplicate_paths!(paths)
    bestPath = find_longest_path(paths)
    w = FASTA.Writer(open("component_$(nbComponent)_longestPath.fasta", "w"))
    best = FASTA.Record("component_$(nbComponent)_longestPath",bestPath.seq)
    write(w, best)
    flush(w)

    w = FASTA.Writer(open("component_$(nbComponent)_allPaths.fasta", "w"))
    for (i,path) in enumerate(paths)
        path = paths[i]
        rec = FASTA.Record("component_$(nbComponent)_path_$(i)",path.seq)
        write(w, rec)
    end
    flush(w)
end
