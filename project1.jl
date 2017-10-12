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
    graphs2search = possibleGraphs(size(data,2))
    #
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
    return Dict(keys[i]=>Names[i] for i = 1:1:l)
end

############################################################
# Function: possibleGraphs(data)
#
# Description: 
############################################################
function possibleGraphs(n::Integer)
    graphs2Search = Set{DiGraph}()
    baseGraph = DiGraph 
    for i = 1:1:n
        
    end
end

#if length(ARGS) != 2
#    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
#end

inputfilename = "small.csv"#ARGS[1]
outputfilename = "short.gph"#ARGS[2]
compute(inputfilename, outputfilename)
