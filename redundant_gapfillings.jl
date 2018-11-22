function merge_gapfillings!(g::MetaBiDiGraph,overlap::Int)
    v = 1
    while v < nv(g)
        merge_redundant_gapfillings!(g,v,"L",overlap)
        merge_redundant_gapfillings!(g,v,"R",overlap)
        v = v+1
    end
    return(g)
end


function merge_redundant_gapfillings!(g::MetaBiDiGraph,startNode::Int,dir::String,overlap::Int)

    if dir =="R" ; rc_strand = "-" ; else rc_strand="+" ; end

    # Starting from a contig node
    if get_prop(g,startNode,:type)=="gapfilling"
        return(g)
    end
    if dir =="R"
        gapfillings = outneighbors(g,startNode)
    else
        gapfillings = inneighbors(g,startNode)
    end

    # Get sequences
    seqs=Dict{Int,String}()
    for node in gapfillings
        if change_dir(g,node,startNode)
            seq = rc(get_prop(g,node,:seq))
        else
            seq = get_prop(g,node,:seq)
        end
        seqs[node] = seq
    end

    # Find unique sequence starts on a 100bp window
    seqStarts = Vector{String}()
    for (node,seq) in seqs
        if length(seq) > 100
            if !(seq[1:100] in seqStarts)
                push!(seqStarts,seq[1:100])
            end
        end
    end
    nbSeq = 0
    for seqStart in seqStarts
        nbSeq += 1
        ref = ""
        breakPos = Dict{Int,Int}()
        refNode=0
        # We search the first diverging positions between 1 seq and all the others
        for node in gapfillings
            seq=seqs[node]
            if length(seq) > 100
                if seq[1:100]==seqStart
                    if length(ref)==0
                        ref=seq
                        refNode = node
                    else
                        breakPos[node] = compare_strings(ref,seqs[node])
                    end
                end
            end
        end
        if length(breakPos)==0
            continue
            #Nothing to compare
        end
        mergePos = minimum(filter(x -> x>100 , collect(values(breakPos)))) # first diverging position
        consensus = ref[1:mergePos-1]
        breakPos[refNode] = 0 # Just to keep track of the node used as ref

        # Create consensus node
        add_vertex!(g)
        if length(seqStarts)==1
            nodeName = get_prop(g,startNode,:name) * "_extended_" * dir
        else
            nodeName = get_prop(g,startNode,:name) * "_extended_" * dir *"(" * string(nbSeq) *")"
        end
        set_prop!(g,nv(g),:name,nodeName)
        set_prop!(g,nv(g),:seq,consensus)
        set_prop!(g,nv(g),:type,"super contig")

        #add_edge!(g,startNode,nv(g))
        if dir =="R"
            add_edge!(g,SimpleBiEdge(startNode,nv(g),"+","+"))
        else
            add_edge!(g,SimpleBiEdge(startNode,nv(g),"-","+"))
        end
        # Link to previous nodes
        for node in keys(breakPos)
            @assert length(list_edges(g,node,startNode))==1
            if change_dir(g,node,startNode)
                if dir == "R"
                    add_edge!(g,SimpleBiEdge(nv(g),node,"+","-"))
                else
                    add_edge!(g,SimpleBiEdge(nv(g),node,"+","+"))
                end
            else
                if dir == "R"
                    add_edge!(g,SimpleBiEdge(nv(g),node,"+","-"))
                else
                    add_edge!(g,SimpleBiEdge(nv(g),node,"+","+"))
                end
            end

            # Shorten previous node
            if !change_dir(g,node,startNode)
                set_prop!(g,node,:seq,get_prop(g,node,:seq)[mergePos+1-overlap:end])
            else
                set_prop!(g,node,:seq,get_prop(g,node,:seq)[1:end-mergePos+overlap])
            end
            # remove old edge
            rem_edge!(g,list_edges(g,node,startNode)[1])
        end
    end
    return(g)
end


function compare_strings(str1,str2)
    # return position of the first difference between str1 and str2
    i=1
    while str1[i]==str2[i] && i < min(length(str1),length(str2))
        i+=1
    end
    return(i)
end



function sliding_score(A, n=5)
  ret = cumsum(A)
  ret[n+1:end] = ret[n+1:end] - ret[1:end-n]
  return ret[n:end-1] / n
end
