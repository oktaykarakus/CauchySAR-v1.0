function [H, PSF] = estimatePSF(RF, psf_size)
% [H, PSF] = estimatePSF(RF, psf_size)
% Estimate psf following the paper: O. V. Michailovich and A. Tannenbaum, “Despeckling of medical
% ultrasound images,?IEEE Transactions on Ultrasonics, Ferroelectrics,
% and Frequency Control, vol. 53, no. 1, pp. 64?8, Jan 2006.
%
% RF is 2D radiofrequency image
% psf_size is [height width] of psf

G = log(abs(fft2(RF)));
deltaG = G - medfilt2(G);
lambda = 0.05*abs(deltaG);
R = sign(deltaG).*max(0, abs(deltaG)-lambda);
GR = G-R;

wname = 'bior3.5';
level = 4;
[C,S] = wavedec2(GR,level,wname);
thr = wthrmngr('dw2ddenoLVL','penalhi',C,S,3);
sorh = 's';
H = exp(wdencmp('lvd',C,S,wname,level,thr,sorh));

if ~isempty(psf_size)
    PSF = abs(otf2psf(H,psf_size));
else
    PSF = abs(otf2psf(H));
end