using BioSequences


function graph2contig(g::MetaDiGraph,outdir::String)
    connectedComponents = weakly_connected_components(g.graph)
    connectedComponents = connectedComponents[sortperm(component_size.(connectedComponents,g),rev=true)]

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
            paths=find_all_paths(g,component[1],"+")
        end
        paths = remove_duplicate_paths!(paths)
        bestPath = find_longest_path(paths)
        w = FASTA.Writer(open(joinpath(outdir,"component_$(nbComponent)_longestPath.fasta"), "w"))
        best = FASTA.Record("component_$(nbComponent)_longestPath",bestPath.seq)
        write(w, best)
        flush(w)

        w = FASTA.Writer(open(joinpath(outdir,"component_$(nbComponent)_allPaths.fasta"), "w"))
        for (i,path) in enumerate(paths)
            path = paths[i]
            rec = FASTA.Record("component_$(nbComponent)_path_$(i)",path.seq)
            write(w, rec)
        end
        flush(w)
    end
end
