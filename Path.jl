import Base.length
import Base.show
import Base.copy

isgapfilling = r".+;.+;len_[0-9]+_qual_[0-9]+_median_cov_[0-9]+"

using MetaGraphs


mutable struct Path
    nodes::Vector{String}
    strands::Vector{String}
    seq::String
    pathName::String

    function Path(g::MetaDiGraph,v::Int,strand::String)  # Constructor of a one node path (initialization)
        p = props(g,v)
        nodes=Vector{String}()
        push!(nodes,p[:name])
        if strand=="-"
            seq = rc(p[:seq])
            pathName = p[:name] * "_Rc"
        else
            seq=p[:seq]
            pathName = p[:name]
        end
        strands=Vector{String}()
        push!(strands,strand)

        new(nodes,strands,seq,pathName)
    end
    function Path(vlist::Vector{String},strands::Vector{String},seq::String,pathName::String) # More general constructor
        new(vlist,strands,seq,pathName)
    end
end

function show(io::IO, p::Path)
    dir = is_directed(g) ? "directed" : "undirected"
    print(io, "Path of $(length(p.nodes)) nodes for a length of $(length(p.seq))")
end

function copy(p::Path)
    Path(deepcopy(p.nodes),deepcopy(p.strands),deepcopy(p.seq),deepcopy(p.pathName))
end

function is_extendable(p,g,kmerSize)
    lastName = last(p.nodes)
    lastNode = collect(filter_vertices(g,:name,lastName))[1]
    if last(p.strands)=="+" ; extendDir="R" ; else ; extendDir ="L" ; end
    nextNodes = neighbors(g,lastNode,extendDir)

    if length(nextNodes) > 1 || length(nextNodes)==0
        return(false)
    else
        nextNode = first(keys(nextNodes))
        dir = first(values(nextNodes))

        # Check that the path is not looping
        if contains(==,p.nodes,get_prop(g,first(keys(nextNodes)),:name))
            return(false)
        end

        # Check that there is no branching
        if dir == "+"
            revNodes = neighbors(g,nextNode,rev_dir(extendDir))
        else
            revNodes = neighbors(g,nextNode,extendDir)
        end
        if length(revNodes) > 1
            return(false)
        else
            extend_path!(p,g,dir,nextNode,kmerSize)
            return(true)
        end
    end
end

function extend_path!(p,g,dir,node,kmerSize)
    seq = get_prop(g,node,:seq)
    pathName = name = get_prop(g,node,:name)
    strand=last(p.strands)
    if dir=="-"
        strand=rev_strand(strand)
    end
    if strand=="-"
        seq=rc(seq)
        pathName = pathName * "_Rc"
    end
    seq=seq[(kmerSize+1):length(seq)]
    push!(p.nodes,name)
    push!(p.strands,strand)
    p.seq = p.seq * seq
    p.pathName = p.pathName * pathName

    return(p)
end


function extend_path!(g::MetaDiGraph,p::Path)
    res = Vector{Path}()
    extended = false
    if p.strands[end]=="+" dir="R" else dir="L" end
    nextNodes = neighbors(g,find_vertex_byname(g,p.nodes[end]),dir)
    for node in keys(nextNodes)
        if get_prop(g,node,:name) in p.nodes # Found a loop
            push!(res,copy(p))
            if get_prop(g,node,:name) == p.nodes[1] # Circular sequence, we should stop there (Is this useful?)
                break
            end
        else
            extended=true
            push!(res,extend_path!(copy(p),g,nextNodes[node],node,kmerSize))
        end
    end
    return(res,extended)
end

function isCyclic(g::MetaDiGraph,p::Path)
    lastNode = find_vertex_byname(g,p.nodes[end])
    firstNode = find_vertex_byname(g,p.nodes[1])

    if p.strands[end]=="-"
        nextNodes=neighbors(g,lastNode,"L")
    else
        nextNodes=neighbors(g,lastNode,"R")
    end
    if firstNode in keys(nextNodes)
        return(true)
    else
        return(false)
    end
end


function find_all_paths(g::MetaDiGraph,node::Int,dir::String)
    paths = [Path(g,node,dir)]
    nbNodes = sum(length.(paths))
    stop=false
    while !stop
        res = Vector{Path}()
        for p in paths
            extendedP = copy(p)
            extendedPaths, extended = extend_path!(g,extendedP)
            if extended
                res = vcat(res,extendedPaths)
            else
                res = vcat(res,p)
            end
        end
        paths = copy(res)
        newNbNodes = sum(length.(paths))
        if newNbNodes == nbNodes
            stop=true
        else
            nbNodes = newNbNodes
        end
    end
    return(paths)
end

function getNames(p::Path)
    return p.nodes
end

function length(p::Path)
    return length(p.nodes)
end

function getNames(v::Vector{Path})
    res = Vector{String}()
    for p in v
        append!(res,getNames(p))
    end
    return(res)
end


function findAllLinearPaths(g::MetaDiGraph,kmerSize::Int)
    v=1
    paths = Vector{Path}()
    while v < nv(g)
        found=false

        ## Find new start node
        name = get_prop(g,v,:name)
        if length(filter(x -> x == name, getNames(paths)))>=1 # path has already be seen
            filter(x -> x==name, getNames(paths))
            v=v+1
            continue
        end

        ## test start node
        if isPathStart(g,v,"R")==true
            p = Path(g,v,"+")
            found = true
        elseif isPathStart(g,v,"L") == true
            p = Path(g,v,"-")
            found = true
        end

        ## Extend path
        if found
            res=true
            while res
                res=is_extendable(p,g,kmerSize)
            end
            push!(paths,p)
        end

        v=v+1
    end
    return(paths)
end

function isPathStart(g::MetaDiGraph,v::Int,dir::String)
    dirNodes = neighbors(g,v,dir)
    oppNodes = neighbors(g,v,rev_dir(dir))

    if length(oppNodes)==1 # It is a valid startpoint if the node before is branching
        prevNode = first(keys(oppNodes))
        if first(values(oppNodes))=="+"
            revNodes = neighbors(g,prevNode,dir)
        else
            revNodes = neighbors(g,prevNode,rev_dir(dir))
        end
        if length(revNodes)==1
            return(false)
        end
    end

    if length(dirNodes)==1
        nextNode = first(keys(dirNodes))
        if first(values(dirNodes))=="+"
            revNodes = neighbors(g,nextNode,rev_dir(dir))
        else
            revNodes = neighbors(g,nextNode,dir)
        end
        if length(revNodes)==1 && first(keys(revNodes))==v
            return(true)
        else
            return(false)
        end
    else
        return(false)
    end
end

function isDeadEnd(g::MetaDiGraph,v::Int)
    if length(neighbors(g,v,"L"))>=1 && length(neighbors(g,v,"R"))==0
        return("-")
    elseif length(neighbors(g,v,"R"))>=1 && length(neighbors(g,v,"L"))==0
        return("+")
    else return(false)
    end
end

# function isCircular(g::MetaDiGraph,v::Int)
#     max_length = 1000
#     dir="R"
#     nextNodes = neighbors(g,v,dir)
#     for (node,orientation) in nextNodes
#         if orientation = "-"
#     end
# end



function merge_path!(g,p)
    add_vertex!(g)
    set_props!(g,nv(g),Dict(:name=>p.pathName, :seq=>p.seq, :type=>"super contig"))


    # Edges from previous nodes
    firstNodeName = p.nodes[1]
    firstNode = collect(filter_vertices(g,:name,firstNodeName))[1]

    if p.strands[1]=="+"
        prevNodes = neighbors(g,firstNode,"L")
    else
        prevNodes = neighbors(g,firstNode,"R")
    end
    for node in keys(prevNodes)
        add_edge!(g,node,nv(g))

        if prevNodes[node]=="+"
            indir = p.strands[1]
        else
            indir = rev_strand(p.strands[1])
        end
        set_prop!(g, Edge(node,nv(g)), :indir, indir)

        set_prop!(g, Edge(node,nv(g)), :outdir, "+")
    end

    # Edges to next nodes
    lastNodeName = last(p.nodes)
    lastNode = collect(filter_vertices(g,:name,lastNodeName))[1]
    if last(p.strands)=="+"
        nextNodes = neighbors(g,lastNode,"R")
    else
        nextNodes = neighbors(g,lastNode,"L")
    end
    for node in keys(nextNodes)
        add_edge!(g,nv(g),node)

        set_prop!(g, Edge(nv(g),node), :indir, "+")
        if nextNodes[node]=="+"
            outdir = last(p.strands)
        else
            outdir = rev_strand(last(p.strands))
        end

        set_prop!(g, Edge(nv(g),node), :outdir, outdir)
    end

    # Remove previous nodes
    g = rem_vertices_byname!(g,getNames(p))

    return(g)
end
function merge_all_linear_paths!(g,LinearPaths)
    for path in LinearPaths
        merge_path!(g,path)
    end
    return(g)
end


function remove_duplicate_paths!(paths::Vector{Path})
    i=1
    while i <= length(paths)
        j=1
        path1=paths[i]
        while j <= length(paths)
            path2=paths[j]
            if path1!=path2
                if path1.nodes == reverse(path2.nodes)
                    filter!(e->e!=path2,paths)
                end
            end
            j+=1
        end
        i+=1
    end
    return(paths)
end

function find_longest_path(paths::Vector{Path})
    maxLength = 0
    local bestPath
    for path in paths
        if length(path.seq) > maxLength
            bestPath = path
            maxLength = length(path.seq)
        end
    end
    return(bestPath)
end
