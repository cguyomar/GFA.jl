
to_compare = collect(keys(neighbors(g,414,"L")))

function merge_redundant_gapfillings!(g::MetaDiGraph,startNode::Int,dir::String)

    overlap=31
    if dir =="R" ; rc_strand = "-" ; else rc_strand="+" ; end

    # Starting from a contig node
    gapfillings = neighbors(g,startNode,dir)

    # Get sequences
    seqs=Dict{Int,String}()
    for node in keys(gapfillings)
        if gapfillings[node]==rc_strand
            seq = rc(get_prop(g,node,:seq))
        else
            seq = get_prop(g,node,:seq)
        end
        seqs[node] = seq
    end

    seqStarts = Vector{String}()
    for node in keys(seqs)
        if length(seqs[node])>100
            if !(seqs[node][1:100] in seqStarts)
                push!(seqStarts,seqs[node][1:100])
            end
        end
    end

    for seqStart in seqStarts
        ref = ""
        breakPos = Dict{Int,Int}()
        print("New \n")
        refNode=0
        for node in keys(seqs)
            seq=seqs[node]
            if length(seq)>100
                if seq[1:100]==seqStart
                    if length(ref)==0
                        ref=seq
                        refNode = node
                        print("Add ref") ; print(refNode)
                        print("\n")
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
        mergePos = minimum(filter(x -> x>100 , collect(values(breakPos))))
        consensus = ref[1:mergePos-1]
        print("\n :");  print(refNode)  ; print("\n")
        breakPos[refNode] = 0 # Just to keep track of the node used as ref

        # Create consensus node
        add_vertex!(g)
        nodeName = get_prop(g,startNode,:name) * "_extended_" * dir
        set_prop!(g,nv(g),:name,nodeName)
        set_prop!(g,nv(g),:seq,consensus)

        add_edge!(g,startNode,nv(g))
        set_prop!(g,startNode,nv(g),:outdir,"+")
        if dir =="R"
            set_prop!(g,startNode,nv(g),:indir,"+")
        else
            set_prop!(g,startNode,nv(g),:indir,"-")
        end

        # Link to old nodes
        for node in keys(breakPos)
            print(node)
            add_edge!(g,nv(g),node)
            set_prop!(g,nv(g),node,:indir,"+")
            if dir == "R"
                set_prop!(g,nv(g),node,:outdir,gapfillings[node])
            else
                set_prop!(g,nv(g),node,:outdir,rev_strand(gapfillings[node]))
            end
            # Shorten old node
            if (gapfillings[node] == "+" && dir =="R") || (gapfillings[node] == "-" && dir =="L")
                set_prop!(g,node,:seq,get_prop(g,node,:seq)[mergePos+1-overlap:end])
            else
                set_prop!(g,node,:seq,get_prop(g,node,:seq)[1:end-mergePos+overlap])
            end
            # remove old edge
            rem_edge!(g,startNode,node)
            rem_edge!(g,node,startNode)
        end
    end
    return(g)
end

# Todo :
# remove too short nodes?



string1 = "ATCGAT"
string2 = "ATCGTA"

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
