function data = rpFast(data, max_recursion_depth, oversampling_ratio, ...
                       decimation_phdata, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2010-2013 Shaun I. Kelly
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% Author: SIK
%
% Email: shaun.kelly@ed.ac.uk
%
% Title: Fast SAR Re-projection
%
% Version: 0.3.0
%
% Description: This function performs the fast re-projection
% operation as described in [1].
%
% Inputs:
% data.Nfft,   Size of the FFT to form the range profile
% data.deltaF, Step size of frequency data (Hz)
% data.minF,   Vector containing the start frequency of each pulse (Hz)
% data.x_mat,  The x-position of each pixel (m)
% data.y_mat,  The y-position of each pixel (m)
% data.z_mat,  The z-position of each pixel (m)
% data.AntX,   The x-position of the sensor at each pulse (m)
% data.AntY,   The y-position of the sensor at each pulse (m)
% data.AntZ,   The z-position of the sensor at each pulse (m)
% data.R0,     The range to scene center (m)
% data.im_final,
%              The complex image value at each pixel
% max_recursion_depth,
%              The number of levels of recursion used
% oversampling_ratio,
%              A ratio of oversampling in range
%              and cross-range, which is a power of two and less than
%              or equal to two to the power of the number of levels
%              of recursion, i.e.
%              oversampling_ratio = 2^N <= 2^max_recursion_depth
% decimation_phdata,
%              decimation in phase history = 1 or
%              decimation in image = 0
% filter_size(optional),
%              Upsampling filter lengths. filter_size(1) specifies
%              the range filter length. filter_size(2) specifies
%              the cross-range filter length. (default value, [41, 41])
% data.beam_center(optional),
%              The Cartesian co-ordinates of the beam center for all
%              or each aperture postion. If it is specified for each
%              postion the dimensions should be 3 * size(data.phdata,
%              2) (default value, [0, 0, 0]) (m, m, m)
% data.antenna_azimuth_beamwidth(optional),
%              The beamwidth of the antenna in azimuth angle.
%              (default value, 2 * pi) (rad)
% data.antenna_polar_beamwidth(optional),
%              The beamwidth of the antenna in polar angle.
%              (default value, 2 * pi) (rad)
%
% Outputs:
% data.phdata, Phase history data (frequency domain) Fast time in
%              rows, slow time in columns
%
% Example:
% data = rpFast(data, max_recursion_depth, oversampling_ratio, ...
%               decimation_phdata, 'filter_size', [41, 41]);
%
% Implementation limitations:
%  * Phase history data must be sampled uniformly in range and
%  cross-range.
%  * Complex image must be defined on approximately square grid.
%
% References:
%
% [1] Kelly, S.I., Rilling, G., Davies, M.E. and Mulgrew B.,
%     "Iterative image formation using fast (re/back)-projection for
%      spotlight-mode SAR,"
%     IEEE Radar Conference, 2011.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global c
global ref_low_pass_filt_1
global ref_low_pass_filt_2
global sub_split_index

c = 299792458;
sub_split_index = [-0.5, -0.5; ...
                   -0.5, +0.5; ...
                   +0.5, -0.5; ...
                   +0.5, +0.5];

% parse input arguments
p = inputParser;
p.addRequired('data', @isstruct);
p.addRequired('max_recursion_depth', @(x) (round(x) == x) && (x >= 0));
p.addRequired('oversampling_ratio', @(x) round(log2(x)) == log2(x) && x > 0);
p.addRequired('decimation_phdata', @(x)  x == 0 || x == 1);
p.addOptional('filter_size', [41, 41], @(x) isvector(x) && length(x) == 2);
p.parse(data, max_recursion_depth, oversampling_ratio, decimation_phdata, ...
        varargin{:});
clear varargin;

c = 299792458;

% set minF elements to the same value
minF = data.minF;
minminF = min(minF);
data.minF = minminF * ones(size(minF));

% make filter odd size
filter_size = 2 * floor(p.Results.filter_size / 2) + 1;

% design reference low-pass filter
[ref_low_pass_filt_1, ref_low_pass_filt_2] = design_ref_filter(filter_size);

data = recursive_rp(data, max_recursion_depth, 0, ...
                    floor(log2(oversampling_ratio)), decimation_phdata);

% restore minF elements and correspondingly correct phdata
r_vec = (c / (2 * data.deltaF * data.Nfft)) * ...
        [-floor(data.Nfft / 2):floor((data.Nfft - 1) / 2)];
for ii = 1:size(data.phdata, 2)
    if abs(minminF - minF(ii)) > eps
        rc = fftshift(ifft(data.phdata(:,ii), data.Nfft));
        rc = rc .* exp(1i * 4 * pi * (minminF - minF(ii)) * r_vec' / c);
        rc = fft(ifftshift(rc));
        data.phdata(:, ii) = rc(1:size(data.phdata, 1));
    end
end
data.minF = minF;

function data = recursive_rp(data, max_recursion_depth, recursion_depth, ...
                             max_oversampling_depth, decimation_phdata)

    global c
    global sub_split_index
    global ref_low_pass_filt_1
    global ref_low_pass_filt_2

    if recursion_depth == max_recursion_depth
        data = rpBasic(data);
%         data = rpBasicAFRLp(data);
    else

        sub_data = data;
        phdata_size = size(data.phdata);
        im_final_size = size(data.im_final);

        recursion_depth = recursion_depth + 1;

        if recursion_depth >= max_oversampling_depth + 1
            % design 2-D low-pass filter
            if decimation_phdata
                low_pass_filt = design_filter(ref_low_pass_filt_1, ...
                                              ref_low_pass_filt_2, ...
                                              phdata_size);
            else
                low_pass_filt = design_filter(ref_low_pass_filt_1, ...
                                              ref_low_pass_filt_2, ...
                                              im_final_size);
            end
        end

        for ii=1:size(sub_split_index, 1)

            if decimation_phdata
                %
                % Decimation in Phase History
                %

                %% image splitting

                % split image into smaller sub-images
                if sub_split_index(ii, 2) > 0
                    sub_im_final_x_range = ...
                        [1:floor(im_final_size(2) / 2)] + ...
                        ceil(im_final_size(2) / 2);
                else
                    sub_im_final_x_range = [1:ceil(im_final_size(2) / 2)];
                end

                if sub_split_index(ii, 1) > 0
                    sub_im_final_y_range = ...
                        [1:floor(im_final_size(1) / 2)] + ...
                        ceil(im_final_size(1) / 2);
                else
                    sub_im_final_y_range = [1:ceil(im_final_size(1) / 2)];
                end

                sub_data.im_final = data.im_final(sub_im_final_y_range, ...
                                                  sub_im_final_x_range);

                sub_data.x_mat = data.x_mat(sub_im_final_y_range, ...
                                            sub_im_final_x_range);
                sub_data.y_mat = data.y_mat(sub_im_final_y_range, ...
                                            sub_im_final_x_range);
                sub_data.z_mat = data.z_mat(sub_im_final_y_range, ...
                                            sub_im_final_x_range);

                % Calculate new scene center
                im_final_center_x_idx = ...
                    ceil((length(sub_im_final_x_range) - 1) / 2) + 1;
                im_final_center_y_idx = ...
                    ceil((length(sub_im_final_y_range) - 1) / 2) + 1;
                im_final_center_x = sub_data.x_mat(1, im_final_center_x_idx);
                im_final_center_y = sub_data.y_mat(im_final_center_y_idx, 1);
                im_final_center_z = sub_data.z_mat(im_final_center_y_idx, ...
                                                   im_final_center_x_idx);

                sub_data.R0 = ((im_final_center_x - data.AntX).^2 + ...
                               (im_final_center_y - data.AntY).^2 + ...
                               (im_final_center_z - data.AntZ).^2).^0.5;
                radial_shift = sub_data.R0 - data.R0;

                if recursion_depth >= max_oversampling_depth + 1
                    %% phase history decimating

                    sub_data.phdata = data.phdata(1:2:end, 1:2:end);

                    sub_data.deltaF = data.deltaF * 2;
                    sub_data.Nfft = ceil(data.Nfft / 2);

                    sub_data.minF = data.minF(1:2:end);
                    sub_data.AntX = data.AntX(1:2:end);
                    sub_data.AntY = data.AntY(1:2:end);
                    sub_data.AntZ = data.AntZ(1:2:end);
                    sub_data.R0 = sub_data.R0(1:2:end);
                end

                %% RECURSIVE CALL
                sub_data = recursive_rp(sub_data, ...
                                        max_recursion_depth, ...
                                        recursion_depth, ...
                                        max_oversampling_depth, ...
                                        decimation_phdata);

                if ii == 1
                    if recursion_depth == max_recursion_depth
                        num_phdata_x_samples = size(sub_data.phdata, 2);
                    end
                    phdata = complex(zeros(phdata_size));
                end

                if recursion_depth >= max_oversampling_depth + 1
                    %% combine deicmated phase histories

                    % up-sample by two

                    % low-pass filter with over-sample by two to remove circular
                    % convolution effects
                    sub_data.phdata = ...
                        ifft2(sub_data.phdata, ...
                              size(low_pass_filt, 1) / 2, ...
                              size(low_pass_filt, 2) / 2);

                    sub_data.phdata = ...
                        [sub_data.phdata, sub_data.phdata; ...
                         sub_data.phdata, sub_data.phdata];

                    sub_data.phdata = fft2(sub_data.phdata .* low_pass_filt);

                    sub_data.phdata = ...
                        sub_data.phdata(1:phdata_size(1), 1:phdata_size(2));
                end

                % focus sub-image phase histories to the new scene
                % reference point by using the fourier shift property
                sub_data.phdata = sub_data.phdata .* ...
                        exp(-1i * 4 * pi * ...
                            (data.minF(1) + ...
                             data.deltaF * [0:(phdata_size(1) - 1)]).' *...
                            radial_shift' / c);

                % combine the sub phdata into the larger phdata
                phdata = phdata + sub_data.phdata;
            else
                %
                % Decimation in Image
                %

                %% phase history splitting

                % split phase history into smaller phase histories
                if sub_split_index(ii, 2) > 0
                    sub_phdata_x_range = ...
                        [1:floor(phdata_size(2) / 2)] + ...
                        ceil(phdata_size(2) / 2);
                else
                    sub_phdata_x_range = [1:ceil(phdata_size(2) / 2)];
                end

                if sub_split_index(ii, 1) > 0
                    sub_phdata_r_range = ...
                        [1:floor(phdata_size(1) / 2)] + ...
                        ceil(phdata_size(1) / 2);
                else
                    sub_phdata_r_range = [1:ceil(phdata_size(1) / 2)];
                end

                sub_data.phdata = ...
                    data.phdata(sub_phdata_r_range, sub_phdata_x_range);

                sub_data.minF = ...
                    data.minF(sub_phdata_x_range) + ...
                    data.deltaF * (sub_phdata_r_range(1) - 1);

                sub_data.Nfft = ceil(data.Nfft / 2);
                sub_data.AntX = data.AntX(sub_phdata_x_range);
                sub_data.AntY = data.AntY(sub_phdata_x_range);
                sub_data.AntZ = data.AntZ(sub_phdata_x_range);
                sub_data.R0 = data.R0(sub_phdata_x_range);

                phdata_center_r_idx = ceil((phdata_size(1) - 1) / 2) + 1;
                phdata_center_x_idx = ceil((phdata_size(2) - 1) / 2) + 1;
                sub_phdata_center_r_idx = ...
                    ceil((length(sub_phdata_r_range) - 1) / 2) + 1;
                sub_phdata_center_x_idx = ...
                    ceil((length(sub_phdata_x_range) - 1) / 2) + 1;

                sub_data.dR0 = ...
                    sqrt((sub_data.AntX(sub_phdata_center_x_idx) - ...
                          data.x_mat).^2 + ...
                         (sub_data.AntY(sub_phdata_center_x_idx) - ...
                          data.y_mat).^2 + ...
                         (sub_data.AntZ(sub_phdata_center_x_idx) - ...
                          data.z_mat).^2) - ...
                    sub_data.R0(sub_phdata_center_x_idx);

                % center image spectrum to new center by using the
                % fourier shift property
                if recursion_depth == 1
                    sub_data.im_final = data.im_final .* ...
                        exp(-1i * 4 * pi * ...
                            (sub_data.minF(1) + ...
                             sub_data.deltaF * ...
                             (sub_phdata_center_r_idx(1) - 1)) * ...
                            sub_data.dR0 / c);
                else
                    if recursion_depth > max_oversampling_depth + 1
                        dR0 = data.dR0(1:2:end, 1:2:end);
                    else
                        dR0 = data.dR0;
                    end

                    sub_data.im_final = data.im_final .* ...
                        exp(-1i * 4 * pi * ...
                            (((sub_data.minF(1) + ...
                               sub_data.deltaF * ...
                               (sub_phdata_center_r_idx(1) - 1)) * ...
                              sub_data.dR0) - ...
                             ((data.minF(1) + ...
                               data.deltaF * ...
                               (phdata_center_r_idx(1) - 1)) * dR0)) / c);
                end

                if recursion_depth >= max_oversampling_depth + 1
                    %% image decimating

                    % low-pass filter with over-sample by two to remove
                    % circular convolution effects
                    sub_data.im_final = ...
                        ifft2(sub_data.im_final, ...
                              size(low_pass_filt, 1), ...
                              size(low_pass_filt, 2));

                    sub_data.im_final = sub_data.im_final .* low_pass_filt;
                    sub_data.im_final = ...
                        sub_data.im_final((1:size(low_pass_filt, 1) / 2), :) + ...
                        sub_data.im_final((size(low_pass_filt, 1) / 2 + 1):end, :);
                    sub_data.im_final = ...
                        sub_data.im_final(:, (1:size(low_pass_filt, 2) / 2)) + ...
                        sub_data.im_final(:, (size(low_pass_filt, 2) / 2 + 1):end);

                    sub_data.im_final = 4 * fft2(sub_data.im_final);

                    sub_data.im_final = ...
                        sub_data.im_final(1:ceil(im_final_size(1) / 2), ...
                                          1:ceil(im_final_size(2) / 2));

                    sub_data.x_mat = data.x_mat(1:2:end, 1:2:end);
                    sub_data.y_mat = data.y_mat(1:2:end, 1:2:end);
                    sub_data.z_mat = data.z_mat(1:2:end, 1:2:end);
                end

                if recursion_depth == max_recursion_depth
                    % remove centering of the sub-image spectrum by using the
                    % fourier shift property
                    if recursion_depth < max_oversampling_depth + 1
                        dR0 = sub_data.dR0;
                    else
                        dR0 = sub_data.dR0(1:2:end, 1:2:end);
                    end

                    sub_data.im_final = sub_data.im_final .* ...
                        exp(1i * 4 * pi * ...
                            (sub_data.minF(1) + ...
                             sub_data.deltaF * ...
                             (sub_phdata_center_r_idx(1) - 1)) * dR0 / c);
                end

                %% RECURSIVE CALL
                sub_data = recursive_rp(sub_data, ...
                                        max_recursion_depth, ...
                                        recursion_depth, ...
                                        max_oversampling_depth, ...
                                        decimation_phdata);

                if ii == 1
                    if recursion_depth == max_recursion_depth
                        num_im_final_x_samples = size(sub_data.im_final, 2);
                        num_im_final_y_samples = size(sub_data.im_final, 1);
                    end
                    phdata = complex(zeros(phdata_size));
                end

                %% combine split phase histories

                % place the sub phdata into the larger phdata
                phdata(sub_phdata_r_range, sub_phdata_x_range) = ...
                    sub_data.phdata;
            end
        end
        data.phdata = phdata;
    end


function [ref_low_pass_filt_1, ...
          ref_low_pass_filt_2] = design_ref_filter(filter_size)

    % range filter
    ref_low_pass_filt_1 = ...
        fir1(filter_size(1) - 1, 0.5, chebwin(filter_size(1)));

    % cross-range filter
    ref_low_pass_filt_2 = ...
        fir1(filter_size(2) - 1, 0.5, chebwin(filter_size(2)));


function low_pass_filt = design_filter(ref_low_pass_filt_1, ...
                                       ref_low_pass_filt_2, ...
                                       data_size)

    % 2-D filter
    filter_size = length(ref_low_pass_filt_1);
    floor_filter_size_2 = floor(filter_size / 2);

    low_pass_filt_r = zeros(2 * floor((data_size(1) + filter_size - 1) / 2), 1);

    low_pass_filt_r(1:(filter_size - floor_filter_size_2)) = ...
        ref_low_pass_filt_1((end - ...
                             (filter_size - floor_filter_size_2) + 1):end);
    low_pass_filt_r((end - floor_filter_size_2 + 1):end) = ...
        ref_low_pass_filt_1(1:floor_filter_size_2);

    filter_size = length(ref_low_pass_filt_2);
    floor_filter_size_2 = floor(filter_size / 2);

    low_pass_filt_x = zeros(2 * floor((data_size(2) + filter_size - 1) / 2), 1);

    low_pass_filt_x(1:(filter_size - floor_filter_size_2)) = ...
        ref_low_pass_filt_2((end - ...
                             (filter_size - floor_filter_size_2) + 1):end);
    low_pass_filt_x((end - floor_filter_size_2 + 1):end) = ...
        ref_low_pass_filt_2(1:floor_filter_size_2);

    low_pass_filt = fft(low_pass_filt_r) * fft(low_pass_filt_x)';


%
% EOF
%
