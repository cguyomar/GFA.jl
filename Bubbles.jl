using MetaGraphs

function pop_bubble!(g,node)
    n1 = neighbors(g,node,"L")
    n2 = neighbors(g,node,"R")
    if length(n1)==length(n2)==1
        if collect(values(n1))[1]=="+"
            l=neighbors(g,collect(keys(n1))[1],"R")
        else
            l=neighbors(g,collect(keys(n1))[1],"L")
        end
        if collect(values(n2))[1]=="+"
            r=neighbors(g,collect(keys(n2))[1],"L")
        else
            r=neighbors(g,collect(keys(n2))[1],"R")
        end

        inter = intersect(r,l)
        if length(inter) > 1
            # get seqs
            seqs = Dict{Int,String}()

            for i in 1:length(inter)
                if !has_edge(g,first(keys(n1)),inter[i][1]) #opposite strand than node
                    seqs[inter[i][1]] = rc(get_prop(g,inter[i][1],:seq))
                else
                    seqs[inter[i][1]] = get_prop(g,inter[i][1],:seq)
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
