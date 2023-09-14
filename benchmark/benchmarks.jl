using BenchmarkTools, Autocorrelations

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
    bdot[i] = @belapsed dotacf($x, $lags)
    bfft[i] = @belapsed fftacf($x, $lags)
    bsmart[i] = @belapsed acf($x, $lags)
end

println("## Scalar timeseries")
display([["N" "dotacf" "fftacf" "acf"]; N bdot bfft bsmart])

#== ## Non-scalar timeseries ==#
for (i,n) in enumerate(N)
    x = [Tuple(randn(2)) for _ in 1:n]
    lags = 0:size(x,1)-1
    bdot[i] = @belapsed dotacf($x, $lags)
    bfft[i] = @belapsed fftacf($x, $lags)
    bsmart[i] = @belapsed acf($x, $lags)
end

println("## Vector timeseries")
display([["N" "dotacf" "fftacf" "acf"]; N bdot bfft bsmart])
