function rem_vertices_byname!(g,nodenames)
    for nodename in nodenames
        to_remove = collect(filter_vertices(g,:name,nodename))[1]
        rem_vertex!(g,to_remove)
    end
    return(g)
end

function find_vertex_byname(g,nodename)
    node = collect(filter_vertices(g,:name,nodename))[1]
    return(node)
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

function compare_nodes(seqs::Dict{Int,String})
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
            if seqs[ref] == seqs[node]
                push!(remove,node)
                foundmatch=true
            else
                aln=pairalign(GlobalAlignment(),seqs[ref],seqs[node],scoremodel)
                if score(aln)/5 > max(length(seqs[ref]),length(seqs[node]))*0.9
                    push!(remove,node)
                    foundmatch=true
                end
            end
        end
        if foundmatch==false
            push!(uniq,node)
        end
    end
    return(remove)
end


function rev_dir(dir::String)
    if dir=="R"
        return("L")
    elseif dir=="L"
        return("R")
    else
        error("format error")
    end
end
function rev_strand(strand::String)
    if strand=="+"
        return("-")
    elseif strand=="-"
        return("+")
    else
        error("format error")
    end
end
