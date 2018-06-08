%% [H, E, beta, gamma] = solveHE(ka, ks, a, s, eBar, uBar, e1, u1, e2, u2, period, lambda, N, M, isTM)
% This program is initially developed by Weijian Yang, and optimized by 
% Vincent Wang and Weijian Yang, in Connie Chang-Hasnain's Group 
% in University of California, Berkeley, 2013
% Paper reference: C. J. Chang-Hasnain and W. Yang, "High contrast gratings for
% integrated optoelectronics," Adv. Opt. Photon. 4, 379-440 (2012).

%{
Copyright (C) Regents of the University of California

This file is part of the high contrast grating solver (HCG Solver)
Written by Weijian Yang, and optimized by Vincent Wang and Weijian Yang
wjyang@eecs.berkeley.edu

HCG Solver is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License 
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}


%%%%%%%%%%% Surface normal light incident
% This function evaluates the matrices H and E, which are the field overlap
% integral between the waveguide array modes inside HCG and the plane waves
% outside the HCG. It essentially describes the properties of the HCG input
% plane and exit plane. H and E are calculated analytically.

% Input parameters (ka and ks are [1 x M] vectors; all others are scalars):
%	ka			transverse wave vectors in the low index air gaps        [1/m]
%	ks			transverse wave vectors in the high index grating bars	 [1/m]
%   period      grating period                                  [m]
%	a			air gap	width									[m]
%	s			grating bar width								[m]
%	eBar		relative permittivity of grating bar			[1]
%	uBar		relative permeability of grating bar			[1]
%	e1			relative permittivity above and below grating	[1]
%	u1			relative permeability above and below grating	[1]
%	e2			relative permittivity between grating bars		[1]
%	u2			relative permeability between grating bars		[1]
%	lambda		wavelength of light	in vacuum					[m]
%	N			number of orders outside grating (total: N)  	[1]
%	M			number of orders inside grating	(total: M)		[1]
%	isTM		1 for TM polarization, 0 for TE					[boolean]

% Output parameters (H and E are [N x M] matrices, 
%                    beta is a [1 x M] vectors,gamma is a [1 x N] vector)
%	H			H field mode overlap integral				    [1]
%	E			E field mode overlap integral			        [1]
%   beta        z-direction mode propagation constants in the waveguide arrays   [1/m]
%   gamma       z-direction mode propagation constants outside the HCG           [1/m]

function [H, E, beta, gamma] = solveHE(ka, ks, a, s, eBar, uBar, e1, u1, e2, u2, period, lambda, N, M, isTM)
	
	n = 1:1:N;			% IMPORTANT: n is an index; _n_ in equations starts from 0, not 1
	m = 1:1:M;
	
	% Precalculate constants
	gamma = conj( 2*pi * sqrt(e1*u1/lambda^2 - (n-1).^2/period^2) );	% (1 x N vector)
	beta = conj( sqrt((2*pi*sqrt(eBar*uBar)/lambda)^2 - ks.^2) );		% (1 x M vector)
	p = (2*pi/period) * (n-1);											% (1 x N vector)
	aa = a/2;
	ss = s/2;
	dif = period - aa;
	
	% NMat(i) and MMat(i) contain all subscript pairs for H and E
	[MMat, NMat] = meshgrid(m, n);
	
	% convert to matrices
	ka = ka(MMat);
	ks = ks(MMat);
	betaMat = beta(MMat);
	gammaMat = gamma(NMat);
	p = p(NMat);
	
	d = ~(NMat-1);
	d = (2-d) / period;
	
	hAir = d .* cos(ks * ss) * 2 .* eval(ka,aa,p,aa);
	hBar = d .* cos(ka * aa) .* (eval(ks,ss,p,dif) - eval(ks,-ss,p,aa));
	
	H = hAir + hBar;
	
	if (isTM)
		E = betaMat ./ gammaMat .* (e1/e2*hAir + e1/eBar*hBar);
	else
		E = gammaMat ./ betaMat .* (u2/u1*hAir + uBar/u1*hBar);
	end
end

% All parameters are N x M matrices or scalars
% 
% sin(KU - PA)     sin(KU + PA)
% ------------  +  ------------ 
%    2(K-P)           2(K+P)
function z = eval(K, U, P, A)
	z = sin(K.*U-P.*A)./(2*(K-P)) + sin(K.*U+P.*A)./(2*(K+P));
	
	assignin('base','K',K);
	assignin('base','U',U);
	assignin('base','P',P);
	assignin('base','A',A);
		
end
