module Autocorrelations

using DSP: conv
using LinearAlgebra: dot
using Statistics: mean

export acf, acf!, fftacf, fftacf!, dotacf, dotacf!

global const THRESHOLD = 1024

"""
    acf(x [, lags]; demean=false, normalize=false)
Evaluate the autocorrelation function of signal `x`.

By default, the acf is evaluated at all available lags `0:size(x,1)-1`,
but arbitrary `lags` can be optionally specified.

**Keywords**
- `demean`: whether to subtract the mean of `x` before evaluating the acf
- `normalize`: whether to normalize the acf to its lag-0 value

When `length(x) < 1024`, a raw computation (`dotacf`) is used,
whereas a FFT-based algorithm is used for larger arrays (`fftacf`).
"""
function acf(x, lags=default_lags(x); demean=false, normalize=false)
    if size(x,1) < THRESHOLD
        dotacf(x, lags; demean=demean, normalize=normalize)
    else
        fftacf(x, lags; demean=demean, normalize=normalize)
    end
end

"""
    acf!(r, x, lags; demean=false, normalize=false)
Evaluate the autocorrelation of `x` at time `lags` in-place and store output in `r`.
See `acf`.
For this in-place version, `lags` *must* be specified, and `length(r)==length(lags)`.
"""
function acf!(r, x, lags; demean=false, normalize=false)
    if size(x,1) < THRESHOLD
        dotacf!(r, x, lags; demean=demean, normalize=normalize)
    else
        fftacf!(r, x, lags; demean=demean, normalize=normalize)
    end
end

#== Autocovariance function with FFT ==#
"""
    fftacf(x [, lags]; demean=false, normalize=false)
Evaluate the autocorrelation of signal `x` with a FFT-based algorithm (see `acf`).
"""
function fftacf(x::AbstractVector{<:Number}, lags=default_lags(x);
    demean::Bool = false, normalize::Bool = false
)
    S = float(eltype(x))
    out = Vector{S}(undef, length(lags))
    fftacf!(out, x, lags; demean=demean, normalize=normalize)
end
function fftacf(x::AbstractVector{T}, lags=default_lags(x);
    demean::Bool = false, normalize::Bool = false
) where {T<:Union{AbstractVector,NTuple}}
    S = float(eltype(eltype(x)))
    out = Vector{S}(undef, length(lags))
    fftacf!(out, x, lags; demean=demean, normalize=normalize)
end

function fftacf!(r::AbstractVector, x::AbstractVector, lags;
    demean::Bool = false, normalize::Bool = false
)
    do_checks(lags, x, r)
    m = length(lags)
    lx = length(x)
    y = demean ? subtract_mean(x) : x
    A = conv(y, reverse(y))
    for k in 1:m
        d = lags[k]
        r[k] = A[d+lx] / (lx-d)
    end
    if normalize
        r0 = A[lx] / lx
        r ./= r0
    end
    return r
end
function fftacf!(r::AbstractVector, x::AbstractVector{T}, lags;
    demean::Bool = false, normalize::Bool = false
) where {T<:Union{AbstractVector,NTuple}}
    do_checks(lags, x, r)
    m = length(lags)
    lx = length(x)
    y = [getindex.(x, i) for i in eachindex(first(x))]
    demean && subtract_mean!(y)
    A = zeros(eltype(r), 2lx-1)
    for s in y
        A .+= conv(s, reverse(s))
    end
    for k in 1:m
        d = lags[k]
        r[k] = A[d+lx] / (lx-d)
    end
    if normalize
        r0 = A[lx] / lx
        r ./= r0
    end
    return r
end

#== Autocovariance function with dot product ==#
"""
    dotacf(x [, lags]; demean=false, normalize=false)
Evaluate the autocorrelation of signal `x` with a dot-product-based
algorithm (see `acf`).
"""
function dotacf(x::AbstractVector{<:Number}, lags=default_lags(x);
    demean::Bool = false, normalize::Bool = false
)
    S = float(eltype(x))
    out = Vector{S}(undef, length(lags))
    dotacf!(out, x, lags; demean=demean, normalize=normalize)
end
function dotacf(x::AbstractVector{T}, lags=default_lags(x);
    demean::Bool = false, normalize::Bool = false
) where {T<:Union{AbstractVector,NTuple}}
    S = float(eltype(eltype(x)))
    out = Vector{S}(undef, length(lags))
    dotacf!(out, x, lags; demean=demean, normalize=normalize)
end

function dotacf!(r::AbstractVector, x::AbstractVector, lags;
    demean::Bool = false, normalize::Bool = false
)
    do_checks(lags, x, r)
    m = length(lags)
    lx = length(x)
    y = demean ? subtract_mean(x) : x
    for k in 1:m
        d = lags[k]
        r[k] = autodot(y, lx, d) / (lx-d)
    end
    if normalize
        r0 = autodot(y, lx, 0) / lx
        r ./= r0
    end
    return r
end


#== Utility functions ==#

autodot(x::AbstractVector{<:Union{Float32,Float64}}, lx::Int, l::Int) = dot(x, 1:(lx-l), x, (1+l):lx)
autodot(x::AbstractVector, lx::Int, l::Int) = dot(view(x, 1:(lx-l)), view(x, (1+l):lx))

function subtract_mean(x::AbstractVector{<:Number})
    S = float(eltype(x))
    out = S.(x)
    subtract_mean!(out)
    return out
end
function subtract_mean(x::AbstractVector)
    S = float(eltype(eltype(x)))
    out = map.(S, x)
    subtract_mean!(out)
    return out
end
function subtract_mean!(x::AbstractVector{<:Number})
    mx = mean(x)
    for k in eachindex(x)
        x[k] -= mx
    end
end
subtract_mean!(x::AbstractVector) = subtract_mean!.(x)

function do_checks(lags, x, r)
    lx = length(x)
    m = length(lags)
    n = length(r)
    if n != m
        error("lags and the output vector must be of the same length")
    end
    if maximum(lags) >= lx
        error("lags must be less than the sample length.")
    end
end

default_lags(x::AbstractVector) = range(0, size(x,1)-1; step=1)

end # module
