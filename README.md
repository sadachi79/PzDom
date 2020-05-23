# extended PzDom model

5/23/2020; Debugged for sorting of labels.

## About exPzDom.jl

The PzDom model, proposed in arXiv:1603.00959v8 [q-bio.PE], exhibits a set of indicators to analyze population dynamics. This model is further developed in https://doi.org/10.1101/780882 and https://doi.org/10.20944/preprints201911.0055.v1. Now we name this extended version of a PzDom model as "exPzDom". In this code of Julia language describing exPzDom, Re(s), "expected sums", Re(v), RRR, E(l) and "threshold" are calculable. 

Re(s) is an indicator related to a covariance of the denoted data, and 0 < Re(s) < 1 means the population is decreasing, while 1 < Re(s) < 2 means increasing, within the category of Gaussian fluctutation. Re(s) = 1 is a neutral situation. Re(s) > 2 means an explosive increase/decrease of population. Further information is available in arXiv:1603.00959v8 [q-bio.PE]. 

"expected sums" is an expected further increase from the particular time point utilized for a calculation. It is caclculable by (∑N)*Im(s)^(-"Number of data groups"), where ∑N is a sum of all the population numbers at a particular time point. Further information is available in arXiv:1603.00959v8 [q-bio.PE]. 

Re(v) is in non-Archimedean sense, larger values mean potentials for the increase in population; smaller values vice versa. Further information is available in https://doi.org/10.1101/780882. Please also refer https://doi.org/10.1371/journal.pone.0179180.

RRR is a calculated potential of a particular phase transition in the particular group.  Further information is available in https://doi.org/10.1101/780882. 

E(l) is a number of external symmetry of an observed group. Further information is available in https://doi.org/10.1101/780882. 

"threshold" is a threshold caclulated by D/lambda_J, where D is a speed of increase for an individual and lambda_J is a Jeans wavelength-like parameter. If the value is smaller than 1, the population would converge, while the value larger than 1 results in divergence. Further information is available in https://doi.org/10.20944/preprints201911.0055.v1. 

## How to use

The input file should be csv file, with the first row being header and the first column being labels for the time points. All the rest should be data, with NaNs replaced by zeros, among non-negative values. At least one data point among a particular time point should be non-zero. Please replace "yourfile.csv" by the path for your csv file, and run the program.  


## Test environment

* macOS High Sierra 10.13.6
* Julia 1.3.1
