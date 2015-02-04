function [registered, output] = efficient_imregister(fixed, moving, subpixel)
% EFFICIENT_IMREGISTER is a convenience wrapper around
% dftregistration.m
%
% Usage: [registered, output] = efficient_imregister(fixed,moving,subpixel)
%
% -- input --
% fixed     : reference image
% moving    : image to be registered
% subpixel  : accuracy scale
% -- output --
% registered : registerted image
% output     : a value vector, error, diffphase, net_row_shift, net_col_shift
%
% March 2014 (c) Michael Hirsch

    
[output Greg] = dftregistration(fft2(fixed), fft2(moving), 1./subpixel);

registered = real(ifft2(Greg));
registered = registered .* (registered>0);