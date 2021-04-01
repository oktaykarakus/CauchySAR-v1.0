%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CauchySAR_v1-SAR Inverse Problem Solving via Cauchy Proximal Splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for inverse problem solution for SAR image super-resolution 
% via the non-convex Cauchy penalty function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some Important Variables
%       ** orgImage             : Original selected input image. It is
%           either HR or LR depending on the selected image.
%
%       ** filterCoef           : Blurring filter coefficients.
%
%       ** dsSize               : Down-sampling ratio. (2: 1/2 downsampling.)
%
%       ** BSNRdb               : Blurred-signal to noise ratio value in
%           decibels. Defines noise strength.
%
%       ** filtered             : Blurred image.
%
%       ** downsampled          : Blurred-and-down-sampled image, i.e. the
%           LR image.
%
%       ** noise_std            : Additive Gaussian noise standard
%           deviation.
%
%       ** monteCarloSize       : The number of Monte Carlo simulations.
%
%       ** maxIter              : Maximum number of iterations of CPS
%           algorithm.
%
%       ** A & AH               : Forward and inverse operators. 
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
%       ** PSNRxxx              : Variables store PSNR values for Noisy and
%           Cauchy reconstructed images. (PSNRbicubic, PSNRnearest and PSNRCauchy)
%
%       ** SSIMxxx              : Variables store SSIM values for Noisy and
%           Cauchy reconstructed images. (SSIMbicubic, SSIMnearest and SSIMCauchy)
%
%       ** RMSExxx              : Variables store RMSE values for Noisy and
%           Cauchy reconstructed images. (RMSEbicubic, RMSEnearest and RMSECauchy)
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

load('CL1_cr.mat');
orgImage = im2double(x2);
filterCoef = fspecial('gaussian', [5 5], 2);
dsSize = 2;                  
BSNRdb = 30;
filtered = imfilter(orgImage, filterCoef, 'circular');
downsampled = downsample(downsample(filtered, dsSize)', dsSize)';
noise_std = norm(downsampled-mean(mean(downsampled)),'fro')/sqrt(size(downsampled, 1)*size(downsampled, 2)*10^(BSNRdb/10));
monteCarloSize = 1;
maxIter = 200;
mu = 1;
gamma = 100*(sqrt(mu)/2);
errorCriterion = 1e-3;
A = @(x) downsample(downsample(imfilter(x, filterCoef, 'circular'), dsSize)', dsSize)';
AH = @(x) imfilter(upsample(upsample(x, dsSize)', dsSize)', filterCoef, 'circular');
for mcCount = 1:monteCarloSize
    degradedImage = downsampled + noise_std*randn(size(downsampled));
    SRImage = imresize(degradedImage, dsSize);
    u  = zeros(size(SRImage));
    oldSRImage = SRImage;
    epsilon = inf;
    for  k = 1:maxIter
        if epsilon(k) > errorCriterion
            u = SRImage- mu * AH(A(SRImage) - degradedImage);
            SRImage = CauchyProx(u, gamma, mu);
        end
        epsilon(k+1) = max(abs( SRImage(:) - oldSRImage(:) )) / max(abs(oldSRImage(:)));
        oldSRImage = SRImage;
    end
    biCubic = imresize(degradedImage, dsSize, 'bicubic');
    Nearest = imresize(degradedImage, dsSize, 'nearest');
    PSNRbicubic(mcCount) = psnr(biCubic(2:end, 2:end), orgImage(2:end, 2:end));
    PSNRnearest(mcCount) = psnr(Nearest(2:end, 2:end), orgImage(2:end, 2:end));
    PSNRCauchy(mcCount) = psnr(SRImage(2:end, 2:end), orgImage(2:end, 2:end));
    
    SSIMbicubic(mcCount) = ssim(biCubic(2:end, 2:end), orgImage(2:end, 2:end));
    SSIMnearest(mcCount) = ssim(Nearest(2:end, 2:end), orgImage(2:end, 2:end));
    SSIMCauchy(mcCount) = ssim(SRImage(2:end, 2:end), orgImage(2:end, 2:end));
    
    RMSEbicubic(mcCount) = sqrt(mean(mean((biCubic(2:end, 2:end) - orgImage(2:end, 2:end)).^2)));
    RMSEnearest(mcCount) = sqrt(mean(mean((Nearest(2:end, 2:end) - orgImage(2:end, 2:end)).^2)));
    RMSECauchy(mcCount) = sqrt(mean(mean((SRImage(2:end, 2:end) - orgImage(2:end, 2:end)).^2)));
    disp(mcCount)
end
figure;
set(gcf, 'Position', [100 100 1400 620])
subplot('Position', [0.0101, 0.25001, 0.23, 0.75])
imshow(orgImage);       
title('Original Image')
subplot('Position', [0.2601, 0.25001, 0.23, 0.75])
imshow(degradedImage);       
title('Degraded Image ')
subplot('Position', [0.5101, 0.25001, 0.23, 0.75])
imshow(biCubic);   
title(['Bicubic Super-resolved (PSNR = ' num2str(mean(PSNRbicubic)) ' dB)'])
subplot('Position', [0.7601, 0.25001, 0.23, 0.75])
imshow(SRImage);  
title(['CPS Super-resolved (PSNR = ' num2str(mean(PSNRCauchy)) ' dB)'])

subplot('Position', [0.0101, 0.01001, 0.23, 0.22])
imshow(orgImage(201:300, 201:450)); 
subplot('Position', [0.2601, 0.01001, 0.23, 0.22])
imshow(degradedImage(101:150, 101:225));   
subplot('Position', [0.5101, 0.01001, 0.23, 0.22])
imshow(biCubic(201:300, 201:450)); 
subplot('Position', [0.7601, 0.01001, 0.23, 0.22])
imshow(SRImage(201:300, 201:450));