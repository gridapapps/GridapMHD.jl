using DrWatson
using DataFrames
using BSON
using CSV

function do_summary(dir)

datadir=projectdir()*"/"*dir
raw = collect_results(datadir)

cols = eachcol(raw)

vals = map(values(cols)) do col
  map(col) do v
    if isa(v,NamedTuple)
      v.max
    elseif isa(v,Tuple)
      first(v)
    else
      v
    end
  end
end

dfmax = DataFrame(vals,keys(cols))
sort!(dfmax,[:np])

rows = eachrow(copy(dfmax))
dict = Dict{Tuple{Int,Int}}{Any}()
for row in rows
  key = (row.nc,row.np)
  if haskey(dict,key)
    prevrow = dict[key]
    for (k,v) in pairs(row)
      # println(k,v)
      if !ismissing(v) && v < prevrow[k]
        prevrow[k] = v
      end
    end
  else
    dict[key] = row
  end
end

df = DataFrame(collect(values(dict)))
sort!(df,order(:np,rev=false))

# mkpath(plotsdir())
# fn = plotsdir("summary.csv")
fn = datadir*"/summary.csv"
CSV.write(fn,df)
df
end

