# Target Function Approximation for Stochastic Circuit Minimization (TFASC)

This project implements the efficient method of target function approximation for stochastic circuit minimization. Given an error bound, the program finds a minimized SC circuit with the approximation error over the target function satisfying the error bound.

Related papers:
- [1]: Exploring Target Function Approximation for Stochastic Circuit Minimization (Chen and Qian, 2020)

## Requirements

- OS: 32-bit Linux (since the tool MVSIS used here cannot be run in 64-bit OS currently)
- gcc
- g++
- make
- libreadline
- ctags
- PERL
- EDA tools: ABC, MVSIS executive files

## Input Format
Content in the "input.txt" file:
```
degree n
precision m
feature vector
```

## Usage



