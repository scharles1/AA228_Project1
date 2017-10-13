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
    i2names = getDict(data)
    graphs2search = possibleGraphs(8)
    for g in graphs2search
       for e in edges(g)
           show(e)
        end
        @printf("\n\n")
    end
    show(length(graphs2search))
end

############################################################
# Function: getDict(data)
#
# Description: 
############################################################
function getDict(data::DataFrame)
    Names = names(data)
    l = length(Names)
    keys = 1:l
    return Dict(keys[i]=>Names[i] for i = 1:l)
end

############################################################
# Function: k2Search(data)
#
# Description: 
############################################################

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

############################################################
# Function: bayesianScore(g1::DiGraph, g2::Digraph)
#
# Description: 
############################################################
function bayesianScore(g1::DiGraph, data::DataFrame) 
    #assume uniform prior with pseudocounts 1
    for n = vertices(g1)
        parents = getParents(g1, n)
        for parent in parents
            nInstances = maximum(data[parent])
            for q = 1:nInstances
                score += lgamma()
            end
        end
    end
end

############################################################
# Function: getParents(g::DiGraph, child::Integer)
#
# Description: 
############################################################
function getParents(g::DiGraph, child::Integer)
    parents = Vector{Integer}()
    for e = edges(g)
        if (dst(e) == child)
           push!(parents, src(e)) 
        end
    end
    return parents
end

############################################################
# Function: getAlpha(data::DataFrame, parents::Vector{Integer})
#
# Description: 
############################################################
function getAlpha(data::DataFrame, i::Integer, j::Integer, k::Integer)
    
end

#if length(ARGS) != 2
#    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
#end

inputfilename = "small.csv"#ARGS[1]
outputfilename = "short.gph"#ARGS[2]
#compute(inputfilename, outputfilename)
g = DiGraph(5)
add_edge!(g,1,2)
add_edge!(g,1,3)
#add_edge!(g,1,4)
add_edge!(g,2,4)
add_edge!(g,3,4)
add_edge!(g,2,3)
data = readtable("small.csv")
numInstances(data,1)