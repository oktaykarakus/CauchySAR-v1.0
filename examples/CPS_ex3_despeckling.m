%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CauchySAR_v1-SAR Inverse Problem Solving via Cauchy Proximal Splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for inverse problem solution for SAR image de-speckling
% via the non-convex Cauchy penalty function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some Important Variables
%       ** imStat               : Variable shows image selection state.
%                                   1 --> Load ready-to-use image
%                                   2 --> Load from the computer
%
%       ** imType               : Variable shows the type of the selected
%           image.
%                                   1 --> SAR - Sea Surface (Speckle-free)
%                                   2 --> SAR - Mountain (Speckle-free)
%                                   3 --> Speckled SAR image
%
%       ** orgImage             : Original selected input image. It is
%           either speckle-free or speckled depending on the selected image.
%
%       ** speckleModel         : Speckle distribution.
%                                   1 --> Gamma distribution.
%                                   2 --> Lognormal distribution.
%                                   3 --> No-noise mimicking.
%
%       ** monteCarloSize       : The number of Monte Carlo simulations.
%
%       ** numLooks             : Number of looks.
%
%       ** noiseSeq             : Generated speckle noise sequence.
%
%       ** speckledImage        : Image after speckle addition.
%
%       ** logSpeckledImage     : Logarithm of speckled image.
%
%       ** waveletMat           : Variable stores wavelet coefficients of
%           the input speckled image.
%
%       ** resultMat            : Variable stores the reconstructed wavelet
%           coefficients via CPS.
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
%       ** reconsWaveCoef       : CPS output.
%
%       ** despeckledImage      : CPS reconstructed de-dspeckled image.
%
%       ** PSNRxxx              : Variables store PSNR values for Noisy and
%           Cauchy reconstructed images. (PSNRNoisy and PSNRCauchy)
%
%       ** SMSExxx              : Variables store S/MSE values for Noisy
%           and Cauchy reconstructed images. (SMSENoisy and SMSECauchy)
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
%% Image selection
imStat = 1;  %                      1 --> Load ready-to-use image
%                                   2 --> Load from the computer
%
imType = 1;  %                      1 --> SAR - Sea Surface (Speckle-free)
%                                   2 --> SAR - Mountain (Speckle-free)
%                                   3 --> Speckled SAR image
%%
if imType == 3
    imStat = 2;
end
if imStat == 1
    if imType == 1
        load('CL1_cr.mat');
    elseif imType == 2
        load('CL5_cr.mat');
    end
    orgImage = double(x2);
elseif imStat == 2
    [file, path] = uigetfile({'*.tif'; '*.bmp'; '*.gif'; '*.jpg'; '*.tiff'}, 'Select file');
    addpath(path);
    orgImage = imread(file);
    if size(orgImage, 3) > 1
        orgImage = mean(orgImage, 3);
    end
end
[r, c] = size(orgImage);
orgImage = orgImage(1:end-mod(r,2), 1:end-mod(c,2));
clear r c

[m, n] = size(orgImage);
speckleModel = 2;
if imType == 3
    speckleModel = 3;
end
monteCarloSize = 1;
numLooks = 5;
maxIter = 200;
mu  = 0.01;
gamma = 10*(sqrt(mu)/2);
errorCriterion = 1e-3;
for mcCount = 1:monteCarloSize
    % calculate the observed image
    if not(imType == 3)
        noiseSeq = generateSpeckle(numLooks, m, n, speckleModel);
        speckledImage = orgImage.*noiseSeq;
    else
        speckledImage = orgImage;
    end
    logSpeckledImage = log(speckledImage);
    waveletMat = cell(1, 4);
    resultMat = cell(1, 4);
    [waveletMat{1, 1},waveletMat{1, 2},waveletMat{1, 3},waveletMat{1, 4}] = dwt2(logSpeckledImage,'db3');
    for layers = 1:4
        waveletCoef = waveletMat{1, layers};
        reconsWaveCoef = waveletCoef;
        u = zeros(size(reconsWaveCoef));
        oldreconsWaveCoef = reconsWaveCoef;
        tic;
        epsilon = inf;
        for  k = 1:maxIter
            if epsilon(k) > 1e-3
                u = reconsWaveCoef - mu*(reconsWaveCoef - waveletCoef);
                reconsWaveCoef = CauchyProx(u, gamma, mu);
            end
            epsilon(k+1) = max(abs( reconsWaveCoef(:) - oldreconsWaveCoef(:) )) / max(abs(oldreconsWaveCoef(:)));
            oldreconsWaveCoef = reconsWaveCoef;
        end
        resultMat{1, layers} = reconsWaveCoef;
    end
    despeckledImage = (exp(idwt2(resultMat{1, 1}, resultMat{1, 2}, resultMat{1, 3}, resultMat{1, 4}, 'db3')));
    PSNRNoisy(mcCount) = psnr(uint8((speckledImage)), uint8(orgImage)); 
    PSNRCauchy(mcCount) = psnr(uint8((despeckledImage)), uint8(orgImage));
    SMSENoisy(mcCount) = 10*log10(sum(orgImage(:).^2)/sum(((speckledImage(:)) - orgImage(:)).^2));
    SMSECauchy(mcCount) = 10*log10(sum(orgImage(:).^2)/sum(((despeckledImage(:)) - orgImage(:)).^2));    
end
figure;
set(gcf, 'Position', [100 100 1400 620])
subplot('Position', [0.0101, 0.25001, 0.3, 0.75])
imshow(uint8(orgImage));       
title('Original Image')
if not(imType == 3)
    subplot('Position', [0.3401, 0.25001, 0.3, 0.75])
    imshow(uint8(speckledImage));      
    title(['Speckled Image (PSNR = ' num2str(mean(PSNRNoisy)) ' dB)'])
end
subplot('Position', [0.6701, 0.25001, 0.3, 0.75])
imshow(uint8(despeckledImage));   
title(['CPS Despeckled (PSNR = ' num2str(mean(PSNRCauchy)) ' dB)'])

subplot('Position', [0.0101, 0.01001, 0.3, 0.22])
imshow(uint8(orgImage(201:300, 201:450)));  
if not(imType == 3)
    subplot('Position', [0.3401, 0.01001, 0.3, 0.22])
    imshow(uint8(speckledImage(201:300, 201:450)));     
end
subplot('Position', [0.6701, 0.01001, 0.3, 0.22])
imshow(uint8(despeckledImage(201:300, 201:450)));  