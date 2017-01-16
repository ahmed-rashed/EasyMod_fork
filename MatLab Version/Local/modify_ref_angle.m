function [theta,theta_ref] = modify_ref_angle(theta_brut)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Angular frame changing with the origin corresponding to the first
%  frequency in the analysed range.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


theta = 2*pi-theta_brut ; % Orientation changing
theta_ref = theta(1,1) ;  % Angle corresponding the the minimum frequency
theta = theta-theta_ref ; % Origin changing
theta = mod(theta,2*pi) ;
