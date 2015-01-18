function [I_noisy, noiseVar] = addnoise(I,SNR,nature)

%//////// Add noise to image
noiseVar = var(nature(:)) / (10^(SNR/10)); % variance of the noise

%##### gaussian
I_noisy = imnoise(I, 'gaussian', 0, noiseVar); % add noise with sigma = v to image I
% equivalent to imnoise(I, 'gaussian', 0, v)
%{
noise = randn(size(I));                 % gaussian distribution noise
noise = noise - mean(noise(:));         % adjust mean to zero
noise = noise./max(abs(vec(noise)));    % scale noise onto [-1,1]
factor = sqrt(v/var(vec(noise)));       % scaler factor which satisfies the desired SNR
noise = noise*factor;                   % 
I_noisy = I + noise;                    % add noise
I_noisy = clip(I_noisy,1,0);            % chop pixel value onto [0,1]
%}

%##### poisson
% I_noisy = imnoise(I, 'poisson');


%//////// rescale
% I_noisy = I_noisy * Imax;
% I_noisy = I_noisy + Imin;

%//////// Show images
% figure
% subplot(1, 2, 1), imshow(I), title('Original image')
% subplot(1, 2, 2), imshow(I_noisy), title('Noisy image, SNR=5db')

end