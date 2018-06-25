function rem_vertices_byname!(g,nodenames)
    for nodename in nodenames
        to_remove = collect(filter_vertices(g,:name,nodename))[1]
        rem_vertex!(g,to_remove)
    end
    return(g)
end


function neighbors(g::MetaDiGraph,node::Int,dir::String)
    if dir == "R"
        dir1="+" ; dir2="-"
    else
        dir1="-" ; dir2="+"
    end

    res=Dict{Int,String}()
    outn = outneighbors(g,node)
    for n in outn
        if get_prop(g,node,n,:indir)==dir1
            if get_prop(g,node,n,:outdir)==get_prop(g,node,n,:indir)
                res[n]="+"
            else
                res[n]="-"
            end
        end
    end

    inn = inneighbors(g,node)
    for n in inn
        if get_prop(g,n,node,:outdir)==dir2
            if get_prop(g,n,node,:outdir)==get_prop(g,n,node,:indir)
                res[n]="+"
            else
                res[n]="-"
            end
        end
    end
    return(res)
end

function compare_nodes(seqs)
    # Compares a set of vertices sequences, and return vertices numbers to delete
    uniq = Int[]
    push!(uniq,first(keys(seqs)))
    remove = Int[]
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1)
    for node in keys(seqs)
        if contains(==,uniq,node)
            continue
        end

        foundmatch=false
        for ref in uniq
            aln=pairalign(GlobalAlignment(),seqs[ref],seqs[node],scoremodel)
            if score(aln)/5 > max(length(seqs[ref]),length(seqs[node]))*0.9
                push!(remove,node)
                foundmatch=true
                continue
            end
        end
        if foundmatch==false
            push!(uniq,node)
        end
    end
    return(remove)
end




function findLinearPath(g::MetaDiGraph,v::Int)
    # Starting from a vertex with no L neighbors
    path = OrientedNode[]
    push!(path,OrientedNode(g,v,"+"))
    nextnodes = neighbors(g,v,"R")
    while (length(nextnodes)==1)

        print("\n")
        if first(nextnodes)[2]=="+"
            push!(path,OrientedNode(g,first(nextnodes)[1],"+"))
            nextnodes = neighbors(g,first(nextnodes)[1],"R")
        elseif first(nextnodes)[2]=="-"
            push!(path,OrientedNode(g,first(nextnodes)[1],"-"))
            nextnodes = neighbors(g,first(nextnodes)[1],"L")
        end
    end
    return(path)
end
