using LightGraphs
using DataFrames

############################################################
# Function: write_gph(dag::DiGraph, idx2names, filename)
#
# Description: Takes a DiGraph, a Dict of index to names 
# and a output filename to write the graph in `gph` format.
############################################################
function write_gph(dag::DiGraph, idx2names::Dict, filename::String)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s, %s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
end

############################################################
# Function: compute(data)
#
# Description:
############################################################
function compute(infile::String, outfile::String)
    data = readtable(infile)
    i2names = get_dict(data)
    graphs2search = possibleGraphs(10)
    for g in graphs2search
       for e in edges(g)
           show(e)
        end
        @printf("\n\n")
    end
    show(length(graphs2search))
end

############################################################
# Function: get_dic(data)
#
# Description: 
############################################################
function get_dict(data::DataFrames.DataFrame)
    Names = names(data)
    l = length(Names)
    keys = 1:l
    return Dict(keys[i]=>Names[i] for i = 1:l)
end

############################################################
# Function: possibleGraphs(data)
#
# Description: 
############################################################
function possibleGraphs(nVertices::Integer)
    graphs2Search = Vector{DiGraph}()
    queue = Vector{DiGraph}()
    unshift!(graphs2Search, DiGraph(nVertices))
    unshift!(queue, DiGraph(nVertices))
    
    while(!isempty(queue))
        curr = pop!(queue)
        if (!hasMarkovEquivalent(curr, graphs2Search))
            unshift!(graphs2Search, copy(curr))
        end
        while(addRandEdge(curr))
           unshift!(queue,copy(curr)) 
        end
    end
    return graphs2Search
end
        
############################################################
# Function: addRandEdge(g::DiGraph)
#
# Description: "Randomly" adds a new directed edge to the
# the graph. Returns true if it was able to. Otherwise,
# it returns false.
############################################################
function addRandEdge(g::DiGraph)
    for i = 1:nv(g)
        for j = 1:nv(g)
            if i == j
                continue
            elseif has_edge(g, Edge(j,i))
                continue
            elseif add_edge!(g,i,j)
                if has_self_loops(g)
                    rem_edge!(g,i,j)
                else
                    return true
                end
            end
        end
    end
    return false
end

############################################################
# Function: hasMarkovEquivalent(g1::DiGraph, v::Vector{Digraph}
#
# Description: 
############################################################
function hasMarkovEquivalent(g1::DiGraph, v::Vector{DiGraph})
    for g2 in v
        if (!markovEquivalent(g1, g2))
            continue
        else
            return true
        end
    end
    return false
end

############################################################
# Function: markovEquivalent(g1::DiGraph, g2::Digraph)
#
# Description: 
############################################################
function markovEquivalent(g1::DiGraph, g2::DiGraph)
    if ne(g1) != ne(g2)
        return false
    end
    for e1 = edges(g1)
       if !has_edge(g2,e1) && !has_edge(g2, reverse(e1))
           return false
       end
       for e2 = edges(g1)
            if e1 == e2
                continue
            elseif dst(e1) == dst(e2)
                if !has_edge(g2,e1) || !has_edge(g2,e2)
                    return false
                end
            end
       end
    end
    return true
end

#if length(ARGS) != 2
#    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
#end

inputfilename = "small.csv"#ARGS[1]
outputfilename = "short.gph"#ARGS[2]
compute(inputfilename, outputfilename)
