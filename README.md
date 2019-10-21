# Iterative-Ensemble-Smoother
An iterative ensemble smoother (iES) based on regularized Levenburg-Marquardt, see the paper "Iterative Ensemble Smoother as an Approximate Solution to a Regularized Minimum-Average-Cost Problem: Theory and Applications ", by Luo et al., SPE-176023-PA, https://doi.org/10.2118/176023-PA

This depository contains an MATLAB implementation of the aforementioned iES, which is mostly used in ensemble-based reservoir data assimilation (also known as history matching) problems. Our main purpose here is to indicate how this iES is actually implemented in our in-house history matching workflow. To this end, we apply this iES to estimate initial conditions of the Lorentzen 96 model, as was done in the paper SPE-176023-PA.  

# A quick start
To run the codes, here I assume that the depository is downloaded into a local computer (and unzipped, if applicable). In this case, please open MATLAB, navigate to the folder "run_example/L96/", and run the following commands

>> [modelInfo,obvInfo,methodInfo] = setupCase();
>> output = iES(modelInfo,obvInfo,methodInfo);

The MATLAB structure variable "output" contains certain information regarding the experiment results. For instance, we can see it by typing

>> output

output = 

       init_obj: 2.8195e+06
    init_objStd: 1.0720e+06
     state_rmse: [1x20 double]
       obs_rmse: [1x20 double]
            obj: 774.7045
         objStd: 9.8481
       exitFlag: [0 0 1]
           iter: 16
    elapsedTime: 5.9519

More experiment results are saved in the sub-folder "run_example/L96/results/".

# An uninformative introduction
The depository contains four folders, namely, "DA", "models", "run_example", "utilities". "DA" contains data assimilation algorithms ("iES.m"); "models" contains models for assimilation ("L96" for the Lorentzen 96 model); "run_example" contains scripts to set up the experiment (e.g., "setupCase.m" under the sub-folder "run_example/L96/"), and is also used to save experiment results (e.g., "run_example/L96/results/"); and "utilities" contains certain scripts needed to run the code.

Based on the template provided in the depository, it would be possible to expand it to other dynamical systems or assimilation algorithms. 

# Disclaimer
This depository is made available and contributed to under the license that include terms that, for the protection of the contributor, make clear that the depository is offered “as-is”, without warranty, and disclaiming liability for damages resulting from using the depository. This guide is no different. The open content license it is offered under includes such terms.

The code may include mistakes, and can’t address every situation. If there is any question, we encourage you to do your own research, discuss with your community or contact me. 

All views expressed here are my own, and do not represent the opinions of any entities with which I have been, am now or will be affiliated with. 
