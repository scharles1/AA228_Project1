using LightGraphs
using DataFrames
using DataStructures
using BayesNets

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
    g1 = DiGraph(8)
    g2 = runSethSearch(g1, data)
    #runK2Search(g1, data)
    show(bayesianScore(g2, data))
    write_gph(g2, i2names, outfile)
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
# Description: Implementing K2 search
############################################################
function runK2Search(g::DiGraph, data::DataFrame)
    bestScore = -10000000000;
    n = nv(g)
    for v = 1:n
        currScore = bayesianScore(g, data)
        prevScore = currScore - 1
        i = 0
        while(currScore > prevScore)
            i += 1
            add_edge!(g, ((v + i) % n) + 1, v)
            if has_self_loops(g)
                rem_edge!(g, ((v + i) % n) + 1, v)
                continue
            end
            prevScore = currScore
            currScore = bayesianScore(g, data)
            if(currScore > bestScore)
                bestScore = currScore
                for edge in edges(g)
                    show(edge)
                    @printf("\n")
                end
            end
            show(currScore)
            @printf("\n\n")
        end
        rem_edge!(g, ((v + i) % n) + 1, v)
    end
end

############################################################
# Function: runSethSearch(data)
#
# Description: Uses a priority queue
############################################################
function runSethSearch(g::DiGraph, data::DataFrame)
    bestScore = -100000000;
    pQ = PriorityQueue(Base.Order.Reverse)
    visitedGraphs = Vector{DiGraph}()
    enqueue!(pQ, copy(g), 0.1)
    n = nv(g)

    while (!isempty(pQ))
        currG = dequeue!(pQ)
        node = rand(1:n)
        for v in vertices(g)
            e = Edge(node, v)
            if has_edge(currG,e) || has_edge(currG,reverse(e)) || (v == node)
                continue
            end
            if !(add_edge!(currG, node, v))
                continue
            end
            if has_self_loops(currG) || hasMarkovEquivalent(currG, visitedGraphs)
                rem_edge!(currG, node, v)
                continue
            end
            currScore = bayesianScore(currG, data)
            @printf("Current Score: %f\n",currScore)
            if(currScore > bestScore)
                bestGraph = copy(currG)
                bestScore = currScore
                show(bestScore)
                @printf("New Best Score: %d\n",bestScore)
                for edge in edges(currG)
                    show(edge)
                    @printf("\n")
                end
                @printf("\n")
            end
            enqueue!(pQ, copy(currG), currScore)
            push!(visitedGraphs, copy(currG))
        end
    end
    return bestGraph
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

############################################################
# Function: bayesianScore(g1::DiGraph, g2::Digraph)
#
# Description: Utilizes the bayesianScoreCompenent and 
#              the getM function
############################################################
function bayesianScore(g::DiGraph, data::DataFrame) 
    #assume uniform prior with pseudocounts 1
    score = 0.0
    for i = vertices(g)
        qTable = getQTable(g, i, data)
        parents = shift!(qTable)
        score += bayesianScoreComponent(data, i, parents, qTable)
    end
    return score
end

############################################################
function bayesianScoreComponent(data::DataFrame, 
                                node::Integer,
                                parents::Vector{Integer},
                                qTable::Vector{Vector{Integer}})
    subscore = 0.0
    nCInstances = maximum(data[node])
    nPInstances = length(qTable)
    if(nPInstances == 0)
        subscore += lgamma(nCInstances)
        subscore -= lgamma(nCInstances + length(data[node]))

        for k = 1:nCInstances
            subscore += lgamma(1 + length(find(x->x == k, data[node])))
        end
    else
        for j = 1:nPInstances
            subscore += lgamma(nCInstances)
            subscore -= lgamma(nCInstances + getM(data, node, j, 0, parents, qTable))
            for k = 1:nCInstances
                subscore += lgamma(1 + getM(data, node, j, k, parents, qTable))
            end
        end
    end
    return subscore
end

############################################################
function getM(data::DataFrame,
              i::Integer,
              j::Integer, 
              k::Integer,
              parents::Vector{Integer},
              qTable::Vector{Vector{Integer}})
    counts = 1:size(data)[1]
    for n = 1:length(parents)
        counts = intersect(counts, find(m->m == qTable[j][n], data[parents[n]]))
    end
    if k == 0
        return length(counts)
    else
        return length(intersect(counts, find(m->m == k, data[i])))
    end
end


############################################################
# Function: getQTable(g::DiGraph, child::Integer)
#
# Description: Extra recursive definition
############################################################
function getQTable(g::DiGraph, child::Integer, data::DataFrame)
    parents = getParents(g, child)
    parentInstances = Vector{Vector{Integer}}()
    for parent in parents
        parentInstance = Vector{Integer}(2)
        parentInstance[1] = parent
        parentInstance[2] = maximum(data[parent])
        push!(parentInstances, parentInstance)
    end
    table = Vector{Vector{Integer}}()
    qi = Vector{Integer}()
    getQTable(table, parentInstances, qi)
    unshift!(table, parents)
    return table
end

############################################################
function getQTable( table::Vector{Vector{Integer}},
                    pInstances::Vector{Vector{Integer}},
                    qi::Vector{Integer})
    if(isempty(pInstances))
        push!(table, copy(qi))
        return
    end
    currInstance = shift!(pInstances)
    for i = 1:currInstance[2]
        push!(qi, i)
        getQTable( table, pInstances, qi)
        pop!(qi)
    end
    unshift!(pInstances, currInstance)
    return
end

############################################################
# Function: getParents(g::DiGraph, child::Integer)
#
# Description: returns the parents of the child in g
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

inputfilename = "small.csv"#ARGS[1]
outputfilename = "small.gph"#ARGS[2]
compute(inputfilename, outputfilename)