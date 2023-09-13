# Autocorrelations.jl

[![Build Status](https://github.com/mastrof/Autocorrelations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mastrof/Autocorrelations.jl/actions/workflows/CI.yml?query=branch%3Amain)

Fast evaluation of autocorrelation functions
```math
f(\tau) = \mathbb{E}\left[X(t)X(t+\tau)\right]
```
for scalar processes $X$ or
```math
f(\tau) = \mathbb{E}\left[\mathbf{X}(t)\cdot\mathbf{X}(t+\tau)\right]
```
for vector processes $\mathbf{X}$.


Specifically, the `acf` function evaluates autocorrelation estimates
of the form
```math
\hat{f}_\tau = \dfrac{1}{T-\tau} \sum_{t=1}^{T-\tau} X_t X_{t+\tau},
\qquad \tau \in [0, T-1]
```
with optional arguments to subtract the mean of the process or to
normalize the resulting autocorrelation function.
