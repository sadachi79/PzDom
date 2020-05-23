using CSV
using DataFrames
using SpecialFunctions
using NLsolve
using Gadfly
using Cairo
using Fontconfig

data = CSV.read("yourfile.csv", header=true, delim=',')

nrow = size(data)[1]
ncol = size(data)[2]

for j in 2:ncol
    data[j] = convert(Array{BigFloat},data[j])
end

@show data

M = []
i = 1
while i <= nrow
    M = append!(M, data[i, 2:end])
    i += 1
end

label = names(data)[2:end]
labels=label
i = 1
while i <= nrow-1
    labels = vcat(labels, label)
    i += 1
end

ld = Any[cat(M[1], labels[1], dims = 1)]
for j in 2:nrow*(ncol-1)
    ld = push!(ld, cat(M[j], labels[j], dims = 1))
    j += 1
end

L = permutedims(reshape(ld, ncol-1, nrow))

sldata = sort!(L, dims=2, by = x -> x[1], rev=true)

sdata = map(x -> x[1], sldata)

slabels = map(x -> x[2], sldata)

rank = zeros(BigFloat, ncol-1, 1)

for k in 1:(ncol-1)
    rank[k] = log(k)
end

Resvalue = zeros(BigFloat, nrow, ncol-1)
for i in 1:nrow
    for j in 1:(ncol-1)
        Resvalue[i,j] = log(sdata[i,1]/sdata[i,j])/rank[j]
    end
end

SM = []
i = 1
while i <= nrow
    M = append!(SM, Resvalue[i, 1:end])
    i += 1
end

SL = []
i = 1
while i <= nrow
    M = append!(SL, slabels[i, 1:end])
    i += 1
end

sld = Any[cat(SM[1], SL[1], dims = 1)]
for j in 2:nrow*(ncol-1)
    sld = push!(sld, cat(SM[j], SL[j], dims = 1))
    j += 1
end

Res0 = permutedims(reshape(sld, ncol-1, nrow))

function nls(func, params...; ini = [0.0])
    if typeof(ini) <: Number
        r = nlsolve((vout,vin)->vout[1]=func(vin[1],params...), [ini])
        v = r.zero[1]
    else
        r = nlsolve((vout,vin)->vout .= func(vin,params...), ini)
        v = r.zero
    end
    return v, r.f_converged
end

function rep(array)
    array = string(array)
    array = replace(array, "(" => "")
    replace(array, ", true)" => "")
end

for i in 1:nrow
    f(x) = zeta(x)-(sum(sdata[i,:])/sdata[i,1])
    Res0[i,1] = cat(parse(BigFloat, rep(nls(f, ini = 10.000001))), slabels[i,1], dims = 1)
end

Res = sort!(Res0, alg=MergeSort, dims=2, by = x -> x[2])

nlRes = map(x -> x[1], Res)

tnlRes = hcat(data[1:end,1], nlRes)

slabels = sort!(slabels, alg=MergeSort, dims=2)

sl = append!([Symbol("Re(s)")], slabels[1,1:end])

Resdf = DataFrame(tnlRes, Symbol.(sl))

@show Resdf

pRes=[]
jp=2
while jp <= ncol
    i=1
    while i <= nrow
        pRes = cat(pRes, cat(cat(Resdf[i,1], Resdf[i,jp], dims=2), sl[jp], dims=2), dims=1)
        i += 1
    end
    jp += 1
end

pl = append!(append!([Symbol("time")], [Symbol("values")]), [Symbol("label")])

pResdf = DataFrame(pRes, Symbol.(pl))

pres = Gadfly.plot(pResdf, x=:"time", y=:"values", color=:"label", Geom.line, Geom.point, Guide.xlabel("time"), Guide.ylabel("Re(s)"))

img = PNG("Re(s)1", 15cm, 10cm)
draw(img, pres)

function df2mm(df::DataFrame)
    n = size(df, 1)
    mm_raw = [fill(1.0, n, 1)]
    mm_name = ["const"]
    for (name, value) in eachcol(df, true)
        if eltype(value) <: Real
            push!(mm_raw, hcat(BigFloat.(value)))
            push!(mm_name, string(name))
        else
            uvalue = unique(value)
            length(uvalue) == 1 && continue
            dvalue = Dict(v=>i for (i, v) in enumerate(uvalue))
            mvalue = zeros(n, length(uvalue))
            for i in 1:n
                mvalue[i, dvalue[value[i]]] = 1.0
            end
            push!(mm_raw, mvalue[:, 2:end])
            append!(mm_name, string.(name, "_", uvalue[2:end]))
        end
    end
    (data=hcat(mm_raw...), names=mm_name)
end

function lm(df, y, xs)
    df = disallowmissing!(dropmissing!(df[:, [y; xs]]))
    yv = BigFloat.(df[!, y])
    xv, xn = df2mm(df[:, [xs;]])
    params = (transpose(xv)*xv)\(transpose(xv)*yv)
    DataFrame(name = xn, estimate=params)
end

for i in 1:nrow
    sdata[i,:] = sort!(sdata[i,:], rev=true)
end

i = 1
lmdf = DataFrame(a = BigFloat[], b = BigFloat[])
while i <= nrow
    df = DataFrame(Rank = rank[:,1], Sdata = sdata[i,:])
    lmdf = push!(lmdf, lm(df, :Sdata, :Rank)[2])
    i += 1
end

logapp = hcat(data[1:end,1], lmdf)

for i in 1:nrow
    logapp[i,3] =  -logapp[i,3]
end

logapp = rename!(logapp, :"x1" => :"log")

@show logapp

mat = [zero(BigFloat) for a in 1:nrow, b in 1:ncol-1]
idata = hcat(data[1:end,1], mat)

Images = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        Images[i,jp] = sum(sdata[i,:])*exp((Resdf[i,jp]*sum(sdata[i,:]))/((ncol-1)*logapp[i,3]))
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("Im(s)")], slabels[1,1:end])

Imsdf = DataFrame(Images, Symbol.(sl))

@show Imsdf

pIms=[]
jp=2
while jp <= ncol
    i=1
    while i <= nrow
        pIms = cat(pIms, cat(cat(Imsdf[i,1], Imsdf[i,jp], dims=2), sl[jp], dims=2), dims=1)
        i += 1
    end
    jp += 1
end

pl = append!(append!([Symbol("time")], [Symbol("values")]), [Symbol("label")])

pImsdf = DataFrame(pIms, Symbol.(pl))

pes = Gadfly.plot(pImsdf, x=:"time", y=:"values", color=:"label", Geom.line, Geom.point, Guide.xlabel("time"), Guide.ylabel("expected sums"))

img = PNG("expected sums1", 15cm, 10cm)
draw(img, pes)

sldata2 = sort!(L, dims=2, by = x -> x[2])

sldata3 = map(x -> x[1], sldata2)

Rev = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        Rev[i,jp] = log(sldata3[i,jp-1])/log(Imsdf[i,jp])
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("Re(v)")], slabels[1,1:end])

Revdf = DataFrame(Rev, Symbol.(sl))

@show Revdf

pRev=[]
jp=2
while jp <= ncol
    i=1
    while i <= nrow
        pRev = cat(pRev, cat(cat(Rev[i,1], Rev[i,jp], dims=2), sl[jp], dims=2), dims=1)
        i += 1
    end
    jp += 1
end

pl = append!(append!([Symbol("time")], [Symbol("values")]), [Symbol("label")])

pRevdf = DataFrame(pRev, Symbol.(pl))

prev = Gadfly.plot(pRevdf, x=:"time", y=:"values", color=:"label", Geom.line, Geom.point, Guide.xlabel("time"), Guide.ylabel("Re(v)"))

img = PNG("Re(v)1", 15cm, 10cm)
draw(img, prev)

Imv = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        Imv[i,jp] = exp((Revdf[i,jp]*sum(sdata[i,:]))/((ncol-1)*logapp[i,3]))
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("Im(v)")], slabels[1,1:end])

Imvdf = DataFrame(Imv, Symbol.(sl))

@show Imvdf

Rev2 = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        Rev2[i,jp] = sldata3[i,jp-1]/Imvdf[i,jp]
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("Re(v2)")], slabels[1,1:end])

Rev2df = DataFrame(Rev2, Symbol.(sl))

@show Rev2df

Imv2 = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        Imv2[i,jp] = exp((Rev2df[i,jp]*sum(sdata[i,:]))/((ncol-1)*logapp[i,3]))
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("Im(v2)")], slabels[1,1:end])

Imv2df = DataFrame(Imv2, Symbol.(sl))

@show Imv2df

pImv2=[]
jp=2
while jp <= ncol
    i=1
    while i <= nrow
        pImv2 = cat(pImv2, cat(cat(Imv2[i,1], Imv2[i,jp], dims=2), sl[jp], dims=2), dims=1)
        i += 1
    end
    jp += 1
end

pl = append!(append!([Symbol("time")], [Symbol("values")]), [Symbol("label")])

pImv2df = DataFrame(pImv2, Symbol.(pl))

pRRR = Gadfly.plot(pImv2df, x=:"time", y=:"values", color=:"label", Geom.line, Geom.point, Guide.xlabel("time"), Guide.ylabel("Im(v2)"))

img = PNG("RRR1", 15cm, 10cm)
draw(img, pRRR)

El = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        El[i,jp] = log(sldata3[i,jp-1])/log(Imvdf[i,jp])
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("E(l)")], slabels[1,1:end])

Eldf = DataFrame(El, Symbol.(sl))

@show Eldf

pEl=[]
jp=2
while jp <= ncol
    i=1
    while i <= nrow
        pEl = cat(pEl, cat(cat(El[i,1], El[i,jp], dims=2), sl[jp], dims=2), dims=1)
        i += 1
    end
    jp += 1
end

pl = append!(append!([Symbol("time")], [Symbol("values")]), [Symbol("label")])

pEldf = DataFrame(pEl, Symbol.(pl))

pel = Gadfly.plot(pEldf, x=:"time", y=:"values", color=:"label", Geom.line, Geom.point, Guide.xlabel("time"), Guide.ylabel("E(l)"))

img = PNG("E(l)1", 15cm, 10cm)
draw(img, pel)

D = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        D[i,jp] = exp(Resdf[i,jp]/logapp[i,3])
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("D")], slabels[1,1:end])

Ddf = DataFrame(D, Symbol.(sl))

@show Ddf

lambda = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        lambda[i,jp] = 2*π*sqrt((logapp[i,3]*D[i,jp])/(sldata3[i,jp-1]*Imsdf[i,jp]))
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("lambda")], slabels[1,1:end])

lambdadf = DataFrame(lambda, Symbol.(sl))

@show lambdadf

threshold = idata

jp=2
while jp <= ncol
    i=1
    while i <= nrow
        threshold[i,jp] = Ddf[i,jp]/lambdadf[i,jp]
        i += 1
    end
    jp += 1
end

sl = append!([Symbol("threshold")], slabels[1,1:end])

thresholddf = DataFrame(threshold, Symbol.(sl))

@show thresholddf

pthreshold=[]
jp=2
while jp <= ncol
    i=1
    while i <= nrow
        pthreshold = cat(pthreshold, cat(cat(threshold[i,1], threshold[i,jp], dims=2), sl[jp], dims=2), dims=1)
        i += 1
    end
    jp += 1
end

pl = append!(append!([Symbol("time")], [Symbol("values")]), [Symbol("label")])

pthresholddf = DataFrame(pthreshold, Symbol.(pl))

pth = Gadfly.plot(pthresholddf, x=:"time", y=:"values", color=:"label", Geom.line, Geom.point, Guide.xlabel("time"), Guide.ylabel("threshold"))

img = PNG("threshold1", 15cm, 10cm)
draw(img, pth)


