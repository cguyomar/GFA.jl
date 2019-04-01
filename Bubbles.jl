function pop_bubble!(g::MetaBiDiGraph,node::Int)
    n1 = inneighbors(g,node)
    n2 = outneighbors(g,node)

    if length(n1)==length(n2)==1 && n1 != n2
        if change_dir(g,node,n1[1])
            l = inneighbors(g,n1[1])
        else
            l = outneighbors(g,n1[1])
        end
        if change_dir(g,node,n2[1]) == true
            r = outneighbors(g,n2[1])
        else
            r = inneighbors(g,n2[1])
        end

        inter = intersect(l,r)
        if length(inter) > 1
            # get seqs
            seqs = Dict{Int,String}()

            for i in 1:length(inter)
                if change_dir(g,node,n1[1]) != change_dir(g,inter[i],n1[1])  #opposite strand than node
                    seqs[inter[i]] = rc(get_prop(g,inter[i],:seq))
                else
                    seqs[inter[i]] = get_prop(g,inter[i],:seq)
                end
            end

            # Compare seqs
            remove = compare_nodes(seqs)
            remove = [get_prop(g,node,:name)  for node in remove]
            g = rem_vertices_byname!(g,remove)
        end
    end
    return(g)
end



function pop_all_bubbles!(g::MetaBiDiGraph)
    v=1
    while v < nv(g)
        print(v)
        nodeName = get_prop(g,v,:name)
        g = pop_bubble!(g,v)
        if nodeName == get_prop(g,v,:name)
            v = v+1
        end
    end
    return(g)
end
