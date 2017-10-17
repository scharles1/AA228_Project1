using LightGraphs
using DataFrames
using DataStructures

############################################################
# Module Variables
############################################################
inputfilename = "large.csv"#ARGS[1]
outputfilename = "large.gph"#ARGS[2]
nNodes = 50
data = readtable(inputfilename)
nInstances = Array([maximum(data[i]) for i = 1:nNodes])
QTable = Dict(i=>Vector{Vector{Integer}}() for i = 1:nNodes)
parents = Dict(i=>Vector{Integer}() for i = 1:nNodes)
data = readtable(inputfilename)
    
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
    i2names = getDict(data)
    g1 = DiGraph(nNodes)    
    bestScore = runRandomK2Search(g1, i2names, -10000000000.0)

    while(bestScore < -3795)
        @printf("Restarting search\n")
        while !remEdge(g1, Edge(rand(1:nv(g1)), rand(1:nv(g1))))
        end
        bestScore = runRandomK2Search(g1, i2names, bestScore)
    end
    #runK2Search(g1, data)
    write_gph(g1, i2names, outfile)
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
# Function: randomk2Search(data)
#
# Description: Implementing K2 search
############################################################
function runRandomK2Search(g::DiGraph, i2names::Dict, bestScore::Float64)
    x = collect(1:nv(g))
    y1 = x[randperm(length(x))]
    y2 = x[randperm(length(x))]
    for v1 in y1
        currScore = bayesianScore(g)
        prevScore = currScore - 1
        for v2 in y2
            if v1 == v2
                continue
            end
            e = Edge(v2,v1)
            if has_edge(g, e) || has_edge(g, reverse(e))
                continue
            end
            if !addEdge(g, e)
                continue
            end
            if is_cyclic(g)
                remEdge(g, e)
                continue
            end
            prevScore = currScore
            currScore = bayesianScore(g)
            @printf("CurrScore: %d\n", currScore)
            if(currScore > bestScore)
                bestScore = currScore
                @printf("New Best Score: %d\n", bestScore)                
                for edge in edges(g)
                    @printf("%s, %s\n", i2names[src(edge)], i2names[dst(edge)])
                end
            end
            @printf("\n")
            if(currScore < prevScore)
                remEdge(g, e)
            end
        end
    end
    return bestScore
end

############################################################
# Function: runSethSearch(data)
#
# Description: Uses a priority queue
############################################################
function runSethSearch(g::DiGraph, data::DataFrame, i2names::Dict)
    bestScore = -100000000;
    pQ = PriorityQueue(Base.Order.Reverse)
    visitedGraphs = Vector{DiGraph}()
    enqueue!(pQ, copy(g), 0.1)
    n = nv(g)
    bestGraph = g
    
    while (!isempty(pQ))
        currG = dequeue!(pQ)
        prevScore = bayesianScore(currG, data)
        node = rand(1:n)
        for v in vertices(g)
            e = Edge(v, node)
            if has_edge(currG,e) || has_edge(currG,reverse(e)) || (v == node)
                continue
            end
            add_edge!(currG, e)
            if is_cyclic(currG) || hasMarkovEquivalent(currG, visitedGraphs)
                rem_edge!(currG, e)
                continue
            end
            currScore = bayesianScore(currG, data)
            @printf("Current Score: %f  PQ_size: %d\n",currScore, length(pQ))
            if(currScore > bestScore)
                bestGraph = copy(currG)
                bestScore = currScore
                @printf("New Best Score: %d\n",bestScore)
                for edge in edges(currG)
                    @printf("%s, %s\n", i2names[src(edge)], i2names[dst(edge)])
                end
                @printf("\n")
                enqueue!(pQ, copy(currG), currScore)
            end
            push!(visitedGraphs, copy(currG))
            if currScore < prevScore
                rem_edge!(currG, e)
            end
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
function bayesianScore(g::DiGraph) 
    #assume uniform prior with pseudocounts 1
    score = 0.0
    for i = vertices(g)
        score += bayesianScoreComponent(i)
    end
    return score
end

############################################################
function bayesianScoreComponent(i::Integer)
    subscore = 0.0
    nCInstances = nInstances[i]
    nPInstances = length(QTable[i])
    if(nPInstances == 0)
        subscore += lgamma(nCInstances)
        subscore -= lgamma(nCInstances + length(data[i]))

        for k = 1:nCInstances
            subscore += lgamma(1 + length(find(x->x == k, data[i])))
        end
    else
        for j = 1:nPInstances
            subscore += lgamma(nCInstances)
            subscore -= lgamma(nCInstances + getM(i, j, 0))
            for k = 1:nCInstances
                subscore += lgamma(1 + getM(i, j, k))
            end
        end
    end
    return subscore
end

############################################################
function getM(i::Integer, j::Integer, k::Integer)
    counts = 1:size(data)[1]
    for n = 1:length(parents[i])
        counts = intersect(counts, find(m->m == QTable[i][j][n], data[parents[i][n]]))
    end
    if k == 0
        return length(counts)
    else
        return length(intersect(counts, find(m->m == k, data[i])))
    end
end

############################################################
# Function: updateQTable(g::DiGraph, child::Integer)
#
# Description: Extra recursive definition
############################################################
function updateQTable(child::Integer)
    global parents
    global nInstances
    global QTable
    ps = parents[child]
    parentInstances = Vector{Vector{Integer}}()
    for p in ps
        parentInstance = Vector{Integer}(2)
        parentInstance[1] = p
        parentInstance[2] = nInstances[p]
        push!(parentInstances, parentInstance)
    end
    table = Vector{Vector{Integer}}()
    qi = Vector{Integer}()
    updateQTable(table, parentInstances, qi)
    QTable[child] = table
end

############################################################
function updateQTable( table::Vector{Vector{Integer}},
                                           pInstances::Vector{Vector{Integer}},
                                           qi::Vector{Integer})
    if(isempty(pInstances))
        push!(table, copy(qi))
        return
    end
    currInstance = shift!(pInstances)
    for i = 1:currInstance[2]
        push!(qi, i)
        updateQTable( table, pInstances, qi)
        pop!(qi)
    end
    unshift!(pInstances, currInstance)
    return
end

############################################################
# Function: addEdge(g::DiGraph)
#
# Description: 
############################################################
function addEdge(g::DiGraph, e::Edge)
    global parents
    if !add_edge!(g, e)
        return false
    end
    p = parents[dst(e)]
    push!(p, src(e))
    updateQTable(dst(e))
    return true
end

############################################################
# Function: remEdge(g::DiGraph)
#
# Description: 
############################################################
function remEdge(g::DiGraph, e::Edge)
    if !rem_edge!(g, e)
        return false
    end
    p = parents[dst(e)]
    deleteat!(p, find(x->x == src(e), p))
    updateQTable(dst(e))
    return true
end

compute(inputfilename, outputfilename)
