

Pre-requisites for PSO:
  -You need "pure" samples for each class (e-, pi-, mu-, etc.).
  -You might need to use data/macro_interaction.sh for extract only events in which the particle interacts with the detector.

Tutorial to use Particle Swarm Optimization (PSO):

1) Clone the repository

2) Configure/Adapt: The python interface manages the communication between the particles and is run on a portal (use screen because of long runtimes).
The Training of the BDTs is done on the batch system and is implemented in Particle.C. This file is also recompiled when you start the PSO now.
    2.1) Remove .pyc files in case you recompile the PSO after any change!
    2.2) This optimization requires Python2.7 & C & C++ & ROOT. The executable SendIFIC250_dEdx_G05_allvars_catA.sh 
    	 is an example of how to run at IFIC, it may need to be adapted.

3) Build a conf file to suit your needs. Check examples in /config/
     play around with the swarm parameters (at least 25 particles recommended)
     choose the Figure of Merit to optimize (at the moment ROC Integral Average)
     fix Kolmogorov-Smirnoff and Anderson-Darling cut values
     specify TMVA factory and method options
     declare coordinate space you want to search
     specify trees and files
     Variables the swarm starts with
     //pool of additional Variables the swarm will try
     3.1) Recommendation: Run one iteration with KS & AD cuts at 0 and then check 
     +If the best configuration have any hyperparameter too close to the limits.
     +Which values of KS & AD ensure safe results (ROC values very close in both Train/Test sets)

4) Start the Optimization with
    python RunPSO.py Example_PSOConfig.txt
    5.1) A better way is using executables, like:
    nohup ./SendIFIC250_dEdx_G05_allvars_catA.sh > nohup_250_dEdx_G05_allvars_catA.log
    
    This way is more stable and straight-forward. Adapt it to yourself!

5) After each iteration the ten best classifiers are writen to PSOResult.txt
   You can also look at it in the Log file
   The best classifier and all necessary information is written to a .conf file

To analysis the results:
1) Get the info for the best value and scan which particle have that performance, do:
       grep -r output "ROCVALUE"
       where "ROCVALUE" is the score of the best performing "particle" (BDT)
2) Copy the TMVA file of that particle into PlotTools/PerformancePlots/ROOTFILES
2.1) Put an appropiate name for it.

3) In PlotTools/PerformancePlots there are the main scripts and macros to study the performance of the classifier.
3.1) One might need to weight the events for a real estimation of the performance: Use data/macro_events.sh for this task
3.2) Use those weights in the macros in PlotTools/PerformancePlots