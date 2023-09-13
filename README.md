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
