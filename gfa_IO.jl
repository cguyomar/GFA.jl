using LightGraphs
using MetaGraphs


function readGFA(infile::String)
    g = MetaDiGraph(PathDiGraph(0))
    lines = readlines(infile)
    for line in lines
        if ismatch(r"S.*",line)
            nodeVal = split(line,"\t")
            add_vertex!(g)
            set_prop!(g, nv(g), :seq, String(nodeVal[3]))
            set_prop!(g, nv(g), :name, String(nodeVal[2]))


        elseif ismatch(r"L.*",line)
            nodeVal = split(line,"\t")
            lv = first(filter_vertices(g,:name,nodeVal[2]))
            rv = first(filter_vertices(g,:name,nodeVal[4]))
            add_edge!(g,lv,rv)
            set_prop!(g,lv,rv,:indir,nodeVal[3])
            set_prop!(g,lv,rv,:outdir,nodeVal[5])
        end
    end
    return(g)
end

function writeToGfa(g::MetaDiGraph,file::String,k::Int)
    open(file, "w") do f

        for vertex in vertices(g)
            p=props(g,vertex)
            write(f,"S\t"*p[:name]*"\t"*p[:seq]*"\n")
        end

        for edge in edges(g)
            p=props(g,edge.src,edge.dst)
            n1=get_prop(g,edge.src,:name)
            n2=get_prop(g,edge.dst,:name)
            write(f,"L\t"*n1*"\t"*p[:indir]*"\t"*n2*"\t"*p[:outdir]*"\t$(k)M\n")
        end
     end
 end
