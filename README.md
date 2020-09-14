# Memetic Genetic Algorithms for Time Series Compression by Piecewise Linear Approximation

This is the accompanying repository to the paper ... published in ....

## Abstract

Time series are sequences of data indexed by time. Such data are collected in various domains, often in massive amounts, such that storing them proves challenging. Thus, time series are commonly stored in a compressed format. An important compression approach is piecewise linear approximation (PLA), which only keeps a small set of time points and interpolates the remainder linearly. Picking a subset of time points such that the PLA minimizes the mean squared error to the original time series is a challenging task, naturally lending itself to heuristics. We propose the piecewise linear approximation genetic algorithm (PLA-GA) for compressing time series by PLA. The PLA-GA is a memetic (μ+λ) GA that makes use of two distinct operators tailored to time series compression. First, we add special individuals to the initial population that are derived using established PLA heuristics. Second, we propose a novel local search operator that greedily improves a compressed time series. We compare the PLA-GA empirically with existing evolutionary approaches and with a deterministic PLA algorithm, known as Bellman's algorithm, that is optimal for the restricted setting of sampling. In both cases, the PLA-GA approximates the original time series better and quicker. Further, it drastically outperforms Bellman's algorithm with increasing instance size with respect to run time until finding a solution of equal or better quality -- we observe speed-up factors between 7 and 100 for instances of 90,000 to 100,000 data points.

This repository contains
- the appendix explaining the operators used in depth
- the measurements
- the source code, which can be used as follows ...

Please note that the timeseries belong to [The UCR Time Series Classification Archive](https://www.cs.ucr.edu/~eamonn/time_series_data_2018/) and to the [UCI Machine Learning Repository](http://archive.ics.uci.edu/ml)
