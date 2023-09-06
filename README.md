# uiNEPv - Locally Unitarily Invariantizable Nonlinear Eigenvector Problems

Dated 		09-06-2023

## Description
This folder contains testing codes and data for the locally unitarily invariantizable NEPv, as used in the paper:

>[1]: *Convergence of SCF for Locally Unitarily Invariantizable NEPv*
by **Ding Lu** and **Ren-Cang Li**, 2022.
Manuscript available at: https://doi.org/10.48550/arXiv.2212.14098


## Contents
- Main functions 
	- *@BuildTrRatio*:	Generate NEPv for the trace ratio optimization problem.
	- *@BuildSumTrRatio*:	Generate NEPv for the sum of trace ratio optimization problem.
	- *@RunSCF*:	Run SCF with level-shifting.
	- *@SolveNEPv*:	Solve the NEPv and check the rate of convergence.

- Examples files
	- *demo_SumTrRatio_X*:	NEPv for the sum of trace ratio optimization. 
	- *demo_TrRatio_X*:	NEPv for the trace ratio optimization.

- Other files:	Supporting functions contained in the private folder.

## Contact 	
For questions, please contact Ding.Lu@uky.edu  
