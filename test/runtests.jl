using Autocorrelations
using Test

@testset "Autocorrelations.jl" begin
    x = 1:10
    f = acf(x)
    @test typeof(f) == Vector{float(eltype(x))}
    @test length(f) == length(x)
    y = 10:-1:1
    u = collect(zip(x,y))
    f = acf(u)
    @test typeof(f) == Vector{float(eltype(eltype(u)))}
    @test length(f) == length(u)
    lags = 1:2:9
    f = acf(x, lags)
    @test length(f) == length(lags)
    lags = 1:2:11
    @test_throws ErrorException acf(x, lags)
    for AcfFunc! in (dotacf!, fftacf!)
        x = 1:10
        lags = 0:size(x,1)-1
        f = acf(x, lags)
        r = Vector{float(eltype(x))}(undef, length(x))
        AcfFunc!(r, x, lags)
        @test r ≈ f
        lags = 0:size(x,1)-2
        @test_throws ErrorException AcfFunc!(r, x, lags)
    end

    @testset "Mathematical correctness" begin
        # compare to analytical result for simple process xₜ=t
        n = 20
        t = 1:n
        x = t
        f = acf(x)
        k = t .- 1
        f_theory = @. (n-k+1) / 2 * (k + (2n-2k+1) / 3)
        @test f ≈ f_theory

        n = 2000
        t = 1:n
        x = t
        f = acf(x)
        k = t .- 1
        f_theory = @. (n-k+1) / 2 * (k + (2n-2k+1) / 3)
        @test f ≈ f_theory

        # compare to raw calculation for process xₜ = sin(ωt)
        n = 100
        t = 1:n
        ω = 2π / 20
        x = @. sin(ω*t)
        f = acf(x)
        k = t .- 1
        f_raw = zeros(n)
        for k in 0:n-1
            for s in 1:(n-k)
                f_raw[1+k] += x[s]*x[s+k]
            end
            f_raw[1+k] /= (n-k)
        end
        @test f ≈ f_raw
    end

    @testset "Equivalence of FFT and dot-based computation" begin
        for x in (1:10, cumsum(randn(500)), sinc.(range(-2π,2π;length=2000)))
            for demean in (true, false), normalize in (true, false)
                a = dotacf(x; demean=demean, normalize=normalize)
                b = fftacf(x; demean=demean, normalize=normalize)
                @test dotacf(x; demean, normalize) ≈ fftacf(x; demean, normalize)
            end
        end
    end
end
