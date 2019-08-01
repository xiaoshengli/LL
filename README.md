# Linear Time Motif Discovery in Time Series (Li-Lin algorithm, LL)

This repository contains the code accompanying the paper, "[Linear Time Motif Discovery in Time Series](https://epubs.siam.org/doi/pdf/10.1137/1.9781611975673.16)" (Xiaosheng Li and Jessica Lin, SDM 2019). This paper presents an algorithm that can find the exact motif of a given time series in a linear expected time complexity. The algorithm is further modified to find all pairs of subsequences whose distances are below a given threshold value. 

The code can also be used to find the closest pair among a set of data points with slight modification. The code is considered suitable for the situation where the data points have high dimension.

## To Compile the Code

Assume using a Linux system:

`g++ -O3 -o LL LL.cpp -std=c++11`

`g++ -O3 -o LL2 LL2.cpp -std=c++11`

## To Run the Code

`./LL [timeSeriesFile] [timeSeriesLength] [motifLength] [g]`

\[timeSeriesFile\] is the name of the time series file to run, the user needs to place the file in the same directory as LL, \[timeSeriesLength\] is the length of the input time series, \[motifLength\] is the motif length, \[g\] is the number of grids. \[g\] supports cascading setting up to 3 values.

`./LL2 [timeSeriesFile] [timeSeriesLength] [motifLength] [threshold] [g]`

\[timeSeriesFile\] is the name of the time series file to run, the user needs to place the file in the same directory as LL2, \[timeSeriesLength\] is the length of the input time series, \[motifLength\] is the motif length, \[threshold\] is the distance threshold for the motif to find, \[g\] is the number of grids.

## Examples

`./LL data.txt 1000000 10 4`

Output:

```
time series length: 1000000 motif length:10
motif distance(Chebyshev): 2.02865e-05
motif: 
630065, 681413
The running time is: 5.31 seconds
```

`./LL data.txt 1000000 10 2 4 6`

Output:

```
time series length: 1000000 motif length:10
motif distance(Chebyshev): 2.02865e-05
motif: 
630065, 681413
The running time is: 12.68 seconds
```

`./LL2 data.txt 1000000 10 0.0001 4`

Output:

```
time series length: 1000000 motif length: 10 threshold: 0.0001 g: 4
number of motif found: 2
motif and distance (Chebyshev): 
330073, 901243, 6.63788e-05
630065, 681413, 2.02865e-05
The running time is: 2.31 seconds
```

## Citation
```
@inproceedings{li2019linear,
  title={Linear Time Motif Discovery in Time Series},
  author={Li, Xiaosheng and Lin, Jessica},
  booktitle={Proceedings of the 2019 SIAM International Conference on Data Mining},
  pages={136--144},
  year={2019},
  organization={SIAM}
}
```
