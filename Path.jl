import Base.length

isgapfilling = r".+;.+;len_[0-9]+_qual_[0-9]+_median_cov_[0-9]+"

using MetaGraphs

mutable struct Path
    nodes::Vector{String}
    strands::Vector{String}
    seq::String
    pathName::String

    function Path(g::MetaDiGraph,v::Int,strand::String)
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
end

function is_extendable(p,g)
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
            extend_path!(p,g,dir,nextNode)
            return(true)
        end
    end
end

function extend_path!(p,g,dir,node)
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
    seq=seq[64:length(seq)] # Should not be hardcoded
    push!(p.nodes,name)
    push!(p.strands,strand)
    p.seq = p.seq * seq
    p.pathName = p.pathName * pathName

    return(p)
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


function findAllLinearPaths(g::MetaDiGraph)
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
                res=is_extendable(p,g)
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
