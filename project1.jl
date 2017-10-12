using LightGraphs
using DataFrames

"""
    write_gph(dag::DiGraph, idx2names, filename)

Takes a DiGraph, a Dict of index to names and a output filename to write the graph in `gph` format.
"""
function write_gph(dag::DiGraph, idx2names, filename)
    open(filename, "w") do io
        for edge in edges(dag)
            @printf(io, "%s, %s\n", idx2names[src(edge)], idx2names[dst(edge)])
        end
    end
end

function compute(infile, outfile)
    data = readtable(infile)
    idx2names = get_dict(data)
    #
    #@printf(io, infile)
end

function get_dic(data)
    names = names(data)
    l = length(names)
    keys = 1:l
    return Dict(keys[i]=>names[i] for i = 1:1:l)
end

if length(ARGS) != 2
    error("usage: julia project1.jl <infile>.csv <outfile>.gph")
end

inputfilename = ARGS[1]
outputfilename = ARGS[2]

compute(inputfilename, outputfilename)
