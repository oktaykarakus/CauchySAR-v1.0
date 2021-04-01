function LoadSARDataRecons(target)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CauchySAR_v1-SAR Inverse Problem Solving via Cauchy Proximal Splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for loading Civilian vehicle data for a specific target. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%       ** target   : Variable representing the target. For this code it
%           can only get 'Camry' value.
%
% Other important variables:
%       ** data     : Structure type variable which includes various 
%           important parameters for the measurements. For details please
%           see below, and also Civilian data set details in [2].
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
%       arXiv preprint arXiv:2005.00657, 2020.
%
% Reference for the data set:
%
% [2] K. E. Dungan, C. Austin, J. Nehrbass, and L. C. Potter, “Civilian vehicle 
%       radar data domes,” inAlgorithms for synthetic aperture radarImagery XVII, 
%       vol. 7699.    International Society for Optics and Photonics, 2010, 
%       p. 76990P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global structSAR
elev = 30;              % What elevation to image (deg)
load([target '_el' num2str(elev) '.0000.mat']);
%% Define input data parameters here
pol = 'VV';             % What polarization to image (HH,HV,VV)
minaz = 100;              % Minimum azimuth angle (degrees)
maxaz = 104;            % Maximum azimuth angle (degrees)

%% Define image parameters here
data.Wx = 10;           % Scene extent x (m)
data.Wy = 10;           % Scene extent y (m)
data.Nfft = 1024;       % Number of samples in FFT
data.Nx = 251;          % Number of samples in x direction
data.Ny = 251;          % Number of samples in y direction
data.x0 = 0;            % Center of image scene in x direction (m)
data.y0 = 0;            % Center of image scene in y direction (m)

%%
% 1 azimuth degree has 16 samples. For 360 degrees, data has 360*16 = 5760 
% samples
idx = 1:5760;
data.hh = data.hh(:, idx);
data.vv = data.vv(:, idx);
data.hv = data.hv(:, idx);
data.azim = data.azim(idx);
newdata.data = data;
% Determine which pulses are located within specified azimuth angles
I = find(and(newdata.data.azim >= minaz, newdata.data.azim <= maxaz));

% Update the phase history for the selected azimuth angles
switch pol
    case{'HH'}
        data.phdata = newdata.data.hh(:,I);
    case{'VV'}
        data.phdata = newdata.data.vv(:,I);
    case{'HV'}
        data.phdata = newdata.data.hv(:,I);
end

% Update other parameters needed for imaging
data.AntAzim = newdata.data.azim(I);
data.AntElev = newdata.data.elev * ones(size(data.AntAzim));
data.freq = newdata.data.FGHz * 1e9;

% Calculate the minimum frequency for each pulse (Hz)
data.minF = min(data.freq)*ones(size(data.AntAzim));

% Calculate the frequency step size (Hz)
data.deltaF = diff(data.freq(1:2));

% Determine the number of pulses and the samples per pulse
[data.K,data.Np] = size(data.phdata);

% Setup imaging grid
data.x_vec = linspace(data.x0 - data.Wx/2, data.x0 + data.Wx/2, data.Nx);
data.y_vec = linspace(data.y0 - data.Wy/2, data.y0 + data.Wy/2, data.Ny);
[data.x_mat,data.y_mat] = meshgrid(data.x_vec,data.y_vec);
data.z_mat = zeros(size(data.x_mat));
data.im_final = zeros(size(data.x_mat));
structSAR = data;