*****************************************************************************************************************
*********************** CauchySAR_v1 - SAR Inverse Problem Solving via Cauchy Proximal Splitting ****************
*****************************************************************************************************************

This package includes the MATLAB source code for using the Cauchy proximal splitting (CPS) algorithm to solve a number of SAR imaging inverse problems including: 

	1) Super-resolution
	3) Image-Formation
	4) De-speckling
	5) Ship Wake Detection.

Specifically, this package includes three folders:

	1) images		: Stores images for the inverse problems examples.
	2) examples		: Stores four scripts for each of the SAR inverse problems.
	
		2.1) CPS_ex1_superresolution.m
		2.2) CPS_ex2_imageFormation.m
		2.3) CPS_ex3_despeckling.m
		2.4) CPS_ex4_shipWakeDetection.m
		
	3) source functions	: Stores the source functions below:
	
		3.1) ConfirmedHalflines.m        
		3.2) norm_visual_SAR.m            
		3.3) estimatePSF.m                
		3.4) radonT.m                     
		3.5) CauchyProx.m                 
		3.6) forwardOp.m                  
		3.7) rpBasicFarField.m            
		3.8) LoadSARDataRecons.m          
		3.9) generateSpeckle.m            
		3.10) rpFast.m                     
		3.11) bpBasicFarField.m            
		3.12) imageNormalize.m             
		3.13) rpFastFarField.m             
		3.14) bpFast.m                     
		3.15) inverseOp.m                  
		3.16) bpFastFarField.m  
		
*****************************************************************************************************************
LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, 
see <https://www.gnu.org/licenses/>.

Copyright (C) Oktay Karakus <o.karakus@bristol.ac.uk> 
		and 
	      Alin Achim <alin.achim@bristol.ac.uk>, 
	      27-03-2021, University of Bristol, UK
*****************************************************************************************************************
REFERENCE

[1] O Karakus, and A Achim. "On Solving SAR Imaging Inverse Problems Using Non-Convex Regularization 
     with a Cauchy-based Penalty"  IEEE Transactions on Geoscience and Remote Sensing, 2020.
arXiv link 	: https://arxiv.org/abs/2005.00657

**** For the CPS Algorithm and its corresponding Matlab package please also refer to:

[2] O Karakus, P Mayo, and A Achim. "Convergence Guarantees for Non-Convex Optimisation with 
     Cauchy-Based Penalties" IEEE Transactions on Signal Processing, 2020.
arXiv link 	: https://arxiv.org/abs/2003.04798

[3] O Karakus, A Achim. (2020): "Cauchy Proximal Splitting (CPS)". 	
https://doi.org/10.5523/bris.15y437loa26cr2nx8gnn3l4hzi 
*****************************************************************************************************************

