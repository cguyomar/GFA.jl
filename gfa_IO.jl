using LightGraphs
using MetaGraphs
using BiDiGraph

isgapfilling = r".+;.+;len_[0-9]+_qual_([0-9]+)_median_cov_[0-9]+"
#gapfillingQual = r".+;.+;len_[0-9]+_qual_(\w[0-9])+_median_cov_[0-9]+"


function readGFA(infile::String)
    isgapfilling = r".+;.+;len_[0-9]+_qual_([0-9]+)_median_cov_[0-9]+"

     file = open(infile,"r")

    g = MetaBiDiGraph(0)
    lines = readlines(infile) # to remove

    while !eof(file)
        line = readline(file)
        if occursin(r"S.*",line)
            nodeVal = split(line,"\t")
            add_vertex!(g)
            set_prop!(g, nv(g), :seq, String(nodeVal[3]))
            set_prop!(g, nv(g), :name, String(nodeVal[2]))
            if occursin(isgapfilling,nodeVal[2])
                set_prop!(g,nv(g),:type,"gapfilling")
                set_prop!(g,nv(g),:qual,parse(Int,match(isgapfilling,nodeVal[2]).captures[1]))
            else
                set_prop!(g,nv(g),:type,"contig")
            end

        elseif occursin(r"L.*",line)
            nodeVal = split(line,"\t")
            lv = first(filter_vertices(g,:name,nodeVal[2]))
            rv = first(filter_vertices(g,:name,nodeVal[4]))
            add_edge!(g,SimpleBiEdge(lv,rv,nodeVal[3],nodeVal[5]))
        end
    end
    return(g)
end

function writeToGfa(g::MetaBiDiGraph,file::String,k::Int)
    open(file, "w") do f

        for vertex in vertices(g)
            p=props(g,vertex)
            write(f,"S\t"*p[:name]*"\t"*p[:seq]*"\n")
        end

        for edge in edges(g)
            p=props(g,edge)
            n1=get_prop(g,edge.src,:name)
            n2=get_prop(g,edge.dst,:name)
            write(f,"L\t"*n1*"\t"*edge.indir*"\t"*n2*"\t"*edge.outdir*"\t$(k)M\n")
        end
     end
 end
