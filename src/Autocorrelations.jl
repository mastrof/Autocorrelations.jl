module Autocorrelations

using DSP: conv
using LinearAlgebra: dot

export acf, fftacf, fftacf!, dotacf, dotacf!

global const THRESHOLD = 1024

acf(x, lags) = size(x,1) < THRESHOLD ? dotacf(x, lags) : fftacf(x, lags)

#== Autocovariance function with FFT ==#
function fftacf(x::AbstractVector{<:Number}, lags=default_lags(x))
    S = float(eltype(x))
    out = Vector{S}(undef, length(lags))
    fftacf!(out, x, lags)
end
function fftacf(x::AbstractVector{T}, lags=default_lags(x)) where {T<:Union{AbstractVector,NTuple}}
    S = float(eltype(eltype(x)))
    out = Vector{S}(undef, length(lags))
    fftacf!(out, x, lags)
end

function fftacf!(r::AbstractVector, x::AbstractVector, lags)
    do_checks(lags, x, r)
    m = length(lags)
    lx = length(x)
    A = conv(x, reverse(x))
    for k in 1:m
        d = lags[k]
        r[k] = A[d+lx]/(lx-d)
    end
    return r
end
function fftacf!(r::AbstractVector, x::AbstractVector{T}, lags) where {T<:Union{AbstractVector,NTuple}}
    do_checks(lags, x, r)
    m = length(lags)
    lx = length(x)
    y = [getindex.(x, i) for i in eachindex(first(x))]
    A = sum([conv(s, reverse(s)) for s in y])
    for k in 1:m
        d = lags[k]
        r[k] = A[d+lx]/(lx-d)
    end
    return r
end

#== Autocovariance function with dot product ==#
function dotacf(x::AbstractVector{<:Number}, lags=default_lags(x))
    S = float(eltype(x))
    out = Vector{S}(undef, length(lags))
    dotacf!(out, x, lags)
end
function dotacf(x::AbstractVector{T}, lags=default_lags(x)) where {T<:Union{AbstractVector,NTuple}}
    S = float(eltype(eltype(x)))
    out = Vector{S}(undef, length(lags))
    dotacf!(out, x, lags)
end

function dotacf!(r::AbstractVector, x::AbstractVector, lags)
    do_checks(lags, x, r)
    m = length(lags)
    lx = length(x)
    for k in 1:m
        d = lags[k]
        r[k] = autodot(x, lx, d) / (lx-d)
    end
    return r
end


#== Utility functions ==#

autodot(x::AbstractVector{<:Union{Float32,Float64}}, lx::Int, l::Int) = dot(x, 1:(lx-l), x, (1+l):lx)
autodot(x::AbstractVector, lx::Int, l::Int) = dot(view(x, 1:(lx-l)), view(x, (1+l):lx))

function do_checks(lags, x, r)
    lx = length(x)
    m = length(lags)
    n = length(r)
    if n != m
        error("lags and the output vector must be of the same length")
    end
    if maximum(lags) < lx
        error("lags must be less than the sample length.")
    end
end

default_lags(x::AbstractVector) = range(0, size(x,1)-1; step=1)

end # module
