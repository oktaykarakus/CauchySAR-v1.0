%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CauchySAR_v1-SAR Inverse Problem Solving via Cauchy Proximal Splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for inverse problem solution for SAR image formation/reconstruction 
% via the non-convex Cauchy penalty function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some Important Variables
%       ** structSAR            : Structure type variable which includes various 
%           important parameters for the measurements. For details please
%           see the function LoadSARDataRecons.m.
%
%       ** structSARMF          : A copy of structSAR for the Matched
%           Filter calculations.
%
%       ** target               : Variable representing the target. For this 
%           code it can only get 'Camry' value.
% 
%       ** phaseHistData        : Input of the CPS algorithm. The phase
%           history measurements.
%       ** A & AH               : Forward and inverse projection operators. 
%
%       ** normdB               : Minimum intensity value in decibels.
%
%       ** mu                   : CPS Step Size
%
%       ** gamma                : The Cauchy penalty function scale
%           parameter.
%
%       ** epsilon              : A vector variable which stores the error 
%           value for each iteration.
%
%       ** errorCriterion       : Target error value.
%
%       ** maxIter              : Maximum number of iterations of CPS
%           algorithm.
%
%       ** CPSRecons            : CPS Resonctructed Image.
%
%       ** normImCauchy         : Normalised CPS reconstructed image.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICENSE
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% Copyright © Oktay Karakus,PhD 
% o.karakus@bristol.ac.uk
% University of Bristol, UK
% July, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE
%
% [1] O Karakus, and A Achim. "On Solving SAR Imaging Inverse Problems Using 
%       Non-Convex Regularization with a Cauchy-based Penalty" 
%       IEEE Transactions on Geoscience and Remote Sensing, 2020.
%
% [2] O Karakus, P Mayo, and A Achim. "Convergence Guarantees for
%       Non-Convex Optimisation with Cauchy-Based Penalties"
%       IEEE Transactions on Signal Processing, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clearvars
close all
clc
cd('..\')
addpath('.\source functions');
addpath('.\examples');
addpath('.\Images');
cd('.\examples')
global structSAR
global structSARMF

target = 'Camry';
LoadSARDataRecons(target)
structSARMF = structSAR;
phaseHistData = structSAR.phdata;
A  = @(x) forwardOp(x);
AH = @(x) inverseOp(x);
normdB = -50;
mu = 1e-8;
gamma = 2*(sqrt(mu)/2);
epsilon = inf;
errorCriterion = 1e-3;
maxIter = 500;
AHy = AH(phaseHistData);
CPSRecons = zeros(size(structSAR.im_final));
iter = 1;
old_CPSRecons = CPSRecons;
structSARMF = bpBasicFarField(structSARMF);
while (iter < maxIter) && (epsilon(iter) > errorCriterion)
    iter = iter + 1;
    grad = AH(A(CPSRecons)) - AHy;
    Z = CPSRecons - mu*grad;
    CPSRecons = CauchyProx(Z, gamma, mu);
    epsilon(iter) = norm(old_CPSRecons(:) - CPSRecons(:)) / norm(CPSRecons(:));
    old_CPSRecons = CPSRecons;
end
normImMF = max(20*log10(abs(structSARMF.im_final)/(max(abs(structSARMF.im_final(:))))), normdB);
normImCauchy = max(20*log10(abs(CPSRecons)/(max(abs(CPSRecons(:))))), normdB);
figure;
set(gcf, 'Position', [100 100 1200 500])
subplot('Position', [0.0501, 0.0901, 0.42, 0.87])
imagesc(structSARMF.x_vec, structSARMF.x_vec, normImMF);
colormap gray
title('MF Reconstruction')
xlabel('x [meters]')
ylabel('y [meters]')
subplot('Position', [0.5101, 0.0901, 0.46, 0.87])
imagesc(structSAR.x_vec, structSAR.x_vec, normImCauchy);
colormap gray
title('CPS Reconstruction')
xlabel('x [meters]')
ylabel('y [meters]')
colorbar