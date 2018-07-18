


function cut_edges!(g,node,dir)
    if dir == "R"
        dir1="+" ; dir2="-"
    else
        dir1="-" ; dir2="+"
    end

    inn = inneighbors(g,node)
    it = 1
    while it < length(inn)
        n = inn[it]
        if get_prop(g,n,node,:outdir)==dir2
            rem_edge!(g,Edge(n,node))
        else
            it+=1
        end
    end

    outn = outneighbors(g,node)
    it = 1
    while it < length(outn)
        n = outn[it]
        if get_prop(g,node,n,:indir)==dir1
            rem_edge!(g,Edge(node,n))
        else
            it+=1
        end
    end
    return(g)
end



function graph_stats(g::MetaDiGraph)
    node_types = countmap(get_type.(vertices(g),g))
    node_types

    contig_length = length(join(get_seq.(collect(filter_vertices(g,:type,"contig")),g)))
    gapfill_length = length(join(get_seq.(collect(filter_vertices(g,:type,"gapfilling")),g)))
    super_length = length(join(get_seq.(collect(filter_vertices(g,:type,"super contig")),g)))


    println("Number of contigs : " * string(node_types["contig"]) * " for a length of " * string(contig_length) * "bp")
    println("Number of gapfilling : " * string(node_types["gapfilling"]) * " for a length of " * string(gapfill_length) * "bp")
    if "super contig" in keys(node_types)
        println("Number of super contigs : " * string(node_types["super contig"]) * " for a length of " * string(gapfill_length) * "bp")
    end

end

function get_type(node,g)
    return(get_prop(g,node,:type))
end
function get_seq(node,g)
    return(get_prop(g,node,:seq))
end

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
                if BioAlignments.score(aln)/5 > max(length(seqs[ref]),length(seqs[node]))*0.9
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
