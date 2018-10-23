using BioSequences
import LightGraphs.ne

function graph2contig(g::MetaDiGraph,outdir::String,overlap::Int)
    connectedComponents = weakly_connected_components(g.graph)
    connectedComponents = connectedComponents[sortperm(component_size.(connectedComponents,g),rev=true)]

    println("Number of connected components : $(length(connectedComponents))")
    nbComponent = 0
    f = open(joinpath(outdir,"component_stats.csv"),"w")
    for (i,component) in enumerate(connectedComponents)
        write(f,string(i)*","*string(length(g,component,overlap))*"\n")

        #found = false
        cycleFound = false
        linearPaths = Vector{Path}()
        cyclePaths = Vector{Path}()
        visited = Set{Int}()

        if length(component)==1
            linearPaths=find_all_paths(g,component[1],"+",kmerSize)
        end

        for node in component
            if isDeadEnd(g,node)==false
                if node in visited
                    continue
                end
                # Search a new cyclic path?
                tmpPaths = find_all_paths(g,node,"+",kmerSize)
                tmpPaths = filter(x -> isCyclic(g,x), tmpPaths)
                if length(tmpPaths)>0
                    cyclePaths = vcat(cyclePaths,tmpPaths)
                    cycleFound = true
                    # add to visited
                    for p in tmpPaths
                        visited = union(visited,find_vertex_byname.(g,p.nodes))
                    end
                end
            else
                if !cycleFound
                    linearPaths = vcat(linearPaths,find_all_paths(g,node,isDeadEnd(g,node)),kmerSize)
                end
            end
        end

        #if found == false # No deadend found -> circle
        #    paths=find_all_paths(g,component[1],"+")
        #end
        if cycleFound
            paths=cyclePaths
        else
            paths=linearPaths
        end
        paths = remove_duplicate_paths!(paths)
        bestPath = find_longest_path(paths)
        w = FASTA.Writer(open(joinpath(outdir,"component_$(i)_longestPath.fasta"), "w"))
        best = FASTA.Record("component_$(i)_longestPath",bestPath.seq)
        write(w, best)
        flush(w)

        w = FASTA.Writer(open(joinpath(outdir,"component_$(i)_allPaths.fasta"), "w"))
        for (j,path) in enumerate(paths)
            path = paths[j]
            rec = FASTA.Record("component_$(i)_path_$(j)",path.seq)
            write(w, rec)
        end
        flush(w)
    end
    close(f)
end


function length(g::MetaDiGraph,component::Array{Int64,1},overlap::Int)
    # Return the total length of a component (by removing overlap length)
    l = 0
    for e in component
        l += length(get_prop(g,e,:seq))
    end
    l = l - ne(g,component)*overlap
    return(l)
end


function ne(g::MetaDiGraph,component::Array{Int64,1})
    # number of edges in a component
    nb=0
    for e in edges(g)
        if e.src in component
            nb += 1
        end
    end
    return(nb)
end
