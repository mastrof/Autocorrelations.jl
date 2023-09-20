using Autocorrelations
using BenchmarkTools
using Printf

open("benchmarks.md", "w") do io
    @printf io "# Autocorrelations.jl benchmarks\n"
    @printf io "(all times in μs)\n\n"
end

## Setup
N = 2 .^ (4:2:16)
b_dummy = @belapsed nothing
bdot = [b_dummy for _ in N]
bfft = [b_dummy for _ in N]
bsmart = [b_dummy for _ in N]

#== ## Scalar timeseries ==#
for (i,n) in enumerate(N)
    x = cumsum(randn(n))
    lags = 0:size(x,1)-1
    bdot[i] = (@belapsed dotacf($x, $lags)) * 1e6
    bfft[i] = (@belapsed fftacf($x, $lags)) * 1e6
    bsmart[i] = (@belapsed acf($x, $lags)) * 1e6
end

open("benchmarks.md", "a") do io
    @printf io "## Scalar timeseries\n"
    @printf io "%-16s %-16s %-16s %-16s\n" "N" "dotacf" "fftacf" "acf"
    for i in eachindex(N)
        @printf io "%-16d %-16.2e %-16.2e %-16.2e\n" N[i] bdot[i] bfft[i] bsmart[i]
    end
    @printf io "\n"
end

#== ## Non-scalar timeseries ==#
for (i,n) in enumerate(N)
    x = [Tuple(randn(2)) for _ in 1:n]
    lags = 0:size(x,1)-1
    bdot[i] = (@belapsed dotacf($x, $lags)) * 1e6
    bfft[i] = (@belapsed fftacf($x, $lags)) * 1e6
    bsmart[i] = (@belapsed acf($x, $lags)) * 1e6
end

open("benchmarks.md", "a") do io
    @printf io "## Non-scalar timeseries\n"
    @printf io "%-16s %-16s %-16s %-16s\n" "N" "dotacf" "fftacf" "acf"
    for i in eachindex(N)
        @printf io "%-16d %-16.2e %-16.2e %-16.2e\n" N[i] bdot[i] bfft[i] bsmart[i]
    end
    @printf io "\n"
end
