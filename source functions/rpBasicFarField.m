function data = rpBasicFarField(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs a basic Backprojection operation.  It is a    %
% different from bpBasic.m as it assumes the data is in the far field. %
% Instead of using the (x,y,z) position of the antenna, this function  %
% uses input azimuth and elevation angles.  The following fields need  %
% to be populated: 
%                                                                      %
% data.Nfft:  Size of the FFT to form the range profile                %
% data.deltaF:  Step size of frequency data (Hz)                       %
% data.minF:  Vector containing the start frequency of each pulse (Hz) %
% data.x_mat:  The x-position of each pixel (m)                        %
% data.y_mat:  The y-position of each pixel (m)                        %
% data.z_mat:  The z-position of each pixel (m)                        %
% data.AntX:  The x-position of the sensor at each pulse (m)           %
% data.AntY:  The y-position of the sensor at each pulse (m)           %
% data.AntZ:  The z-position of the sensor at each pulse (m)           %
% data.R0:  The range to scene center (m)                              %
% data.phdata:  Phase history data (frequency domain)                  %
%               Fast time in rows, slow time in columns                %
%                                                                      %
% The output is:                                                       %
% data.im_final:  The complex image value at each pixel                %
%                                                                      %
% Written by LeRoy Gorham, Air Force Research Laboratory, WPAFB, OH    %
% Email:  leroy.gorham@wpafb.af.mil                                    %
% Date Released:  8 Apr 2010                                           %
%                                                                      %
% Gorham, L.A. and Moore, L.J., "SAR image formation toolbox for       %
%   MATLAB,"  Algorithms for Synthetic Aperture Radar Imagery XVII     %
%   7669, SPIE (2010).                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define speed of light (m/s)
c = 299792458;

% Determine the size of the phase history data
data.K = size(data.phdata,1);  % The number of frequency bins per pulse
data.Np = size(data.phdata,2); % The number of pulses

% Determine the azimuth angles of the image pulses (radians)
data.AntAz = sort(data.AntAzim*pi/180);

% Determine the average azimuth angle step size (radians)
data.deltaAz = abs(mean(diff(data.AntAz)));

% Determine the total azimuth angle of the aperture (radians)
data.totalAz = max(data.AntAz) - min(data.AntAz);

% Determine the maximum wavelength (m)
data.maxLambda = c / (mean(data.minF) + data.deltaF * data.K);

% Determine the maximum scene size of the image (m)
data.maxWr = c/(2*data.deltaF);   
data.maxWx = data.maxLambda/(2*data.deltaAz);

% Determine the resolution of the image (m)
data.dr = c/(2*data.deltaF*data.K);
data.dx = data.maxLambda/(2*data.totalAz);

% Display maximum scene size and resolution
% fprintf('Maximum Scene Size:  %.2f m range, %.2f m cross-range\n',data.maxWr,data.maxWx);
% fprintf('Resolution:  %.2fm range, %.2f m cross-range\n',data.dr,data.dx);

% Calculate the range to every bin in the range profile (m)
data.r_vec = linspace(-data.Nfft/2,data.Nfft/2-1,data.Nfft)*data.maxWr/data.Nfft;

% Initialize the phase histories with all zero values
data.phdata = zeros(size(data.phdata));

% Set up a vector to keep execution times for each pulse (sec)
t = zeros(1,data.Np);

% Loop through every pulse
for ii = 1:data.Np

    % Display status of the imaging process
%     if ii > 1 && mod(ii,100)==0
%         t_sofar = sum(t(1:(ii-1)));
%         t_est = (t_sofar*data.Np/(ii-1)-t_sofar)/60;
%         fprintf('Pulse %d of %d, %.02f minutes remaining\n',ii,data.Np,t_est);
%     else
%         fprintf('Pulse %d of %d\n',ii,data.Np);
%     end
%     tic

    % Calculate differential range for each pixel in the image (m)
    dR = data.x_mat * cosd(data.AntElev(ii)) * cosd(data.AntAzim(ii)) + ...
        data.y_mat * cosd(data.AntElev(ii)) * sind(data.AntAzim(ii)) + ...
        data.z_mat * sind(data.AntElev(ii));

    % Calculate phase correction for image
    phCorr = exp(-1i*4*pi*data.minF(ii)/c*dR);

    % Determine which pixels fall within the range swath
    I = find(and(dR > min(data.r_vec), dR < max(data.r_vec)));

    % Update the phase histories using linear interpolation
    rc = my_interp(dR(I), data.im_final(I) .* phCorr(I), data.r_vec);
    rc = fft(fftshift(rc),data.Nfft); % fft(fftshift(rc));
    %%% Correction
    temp = fftshift(rc);
%     rc = temp(floor((data.Nfft-data.K)/2)+1:3*floor((data.Nfft-data.K)/2));
    rc = temp(floor((data.Nfft-data.K)/2)+1:floor((data.Nfft+data.K)/2));
    %%%
%     data.phdata(:, ii) = rc(1:size(data.phdata, 1));
    data.phdata(:, ii) = rc;

    % Determine the execution time for this pulse
%     t(ii) = toc;
end

return

function rc = my_interp(dR, im, r_vec)

rc = zeros(size(r_vec))';

r_diff_recip = 1 / (r_vec(2) - r_vec(1));
r_idx = (dR(:) - r_vec(1)) * r_diff_recip + 1;
r_idx_minus = floor(r_idx);
r_idx_plus = r_idx_minus + 1;

im_minus = im(:) .* (r_idx_plus - r_idx);
im_plus = im(:) .* (r_idx - r_idx_minus);

for i=1:length(dR)
    rc(r_idx_minus(i)) = rc(r_idx_minus(i)) + im_minus(i);
    rc(r_idx_plus(i)) = rc(r_idx_plus(i)) + im_plus(i);
end
