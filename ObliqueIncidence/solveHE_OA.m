%% function [H, E, beta, gamma] = solveHE_OA(ka, ks, a, s, period, eBar, uBar, e1, u1, e2, u2, lambda, theta, N, M, isTM)
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


% This function evaluates the matrices H and E, which are the field overlap
% integral between the waveguide array modes inside HCG and the plane waves
% outside the HCG. It essentially describes the properties of the HCG input
% plane and exit plane. H and E are calculated analytically.

% Input parameters (ka and ks are [1 x M] vectors; all others are scalars):
%	ka			transverse wave vectors in the low index air gaps        [1/m]
%	ks			transverse wave vectors in the high index grating bars	 [1/m]
%	a			air gap	width									[m]
%	s			grating bar width								[m]
%	eBar		relative permittivity of grating bar			[1]
%	uBar		relative permeability of grating bar			[1]
%	e1			relative permittivity above and below grating	[1]
%	u1			relative permeability above and below grating	[1]
%	e2			relative permittivity between grating bars		[1]
%	u2			relative permeability between grating bars		[1]
%	lambda		wavelength of light	in vacuum					[m]
%	theta		incident light angle							[rad]
%	N			number of orders outside grating (total: 2N+1)	[1]
%	M			number of orders inside grating	(total: M)		[1]
%	isTM		1 for TM polarization, 0 for TE					[boolean]

% Output parameters (H and E are [(2N+1) x M] matrices, 
%                    beta is a [1 x M] vectors,gamma is a [1 x (2N+1)] vector)
%	H			H field mode overlap integral				    [1]
%	E			E field mode overlap integral			        [1]
%   beta        z-direction mode propagation constants in the waveguide arrays   [1/m]
%   gamma       z-direction mode propagation constants outside the HCG           [1/m]

function [H, E, beta, gamma] = solveHE_OA(ka, ks, a, s, period, eBar, uBar, e1, u1, e2, u2, lambda, theta, N, M, isTM)
	n1 = sqrt(e1*u1);
	k0 = 2*pi/lambda;
	kx0 = n1*k0 * sin(theta);
	
	% define constants
	if (isTM)
		NN = e2 / eBar;
	else
		NN = 1;
	end
	EE = exp(1i*kx0*period);
	
	sa = sin(ka * a/2);
	ss = sin(ks * s/2);
	ca = cos(ka * a/2);
	cs = cos(ks * s/2);
	ss2 = sin(ks * s);
	cs2 = cos(ks * s);
	
	A = NN*ks.*sa + EE*ka.*ca.*ss2 + EE*NN*ks.*sa.*cs2;
	B = NN*ks.*ca + EE*ka.*sa.*ss2 - EE*NN*ks.*ca.*cs2;
	C = ka.*ss.*(EE+2*ca.^2-1) + 2*NN*ks.*ca.*sa.*cs;
	D = -ka.*cs.*(EE-2*ca.^2+1) - 2*NN*ks.*ca.*sa.*ss;
	
	n = -N:1:N;
	kxn = kx0 + 2*pi*n/period;
	gamma = conj(sqrt(e1*u1*k0^2 - kxn.^2));
	beta = conj(sqrt(eBar*uBar*k0^2-ks.^2));
	
	[mMat, nMat] = meshgrid(1:M, 1:2*N+1);
	A = A(mMat);
	B = B(mMat);
	C = C(mMat);
	D = D(mMat);
	kxn = kxn(nMat);
	gammaMat = gamma(nMat);
	betaMat = beta(mMat);
	ka = ka(mMat);
	ks = ks(mMat);
	ca = ca(mMat);
	cs = cs(mMat);
	sa = sa(mMat);
	ss = ss(mMat);
	
	divAir = period * (ka.^2 - kxn.^2);
	divBar = period * (ks.^2 - kxn.^2);

	e0 = 1;
	ea = exp(1i*kxn*a);
	ep = exp(1i*kxn*period);
	hAir = eval(kxn, ka, A, B, e0, ea, ca, sa) ./ divAir;
	hBar = eval(kxn, ks, C, D, ea, ep, cs, ss) ./ divBar;	
	
	H = hAir + hBar;

	if (isTM)
		E = betaMat ./ gammaMat .* (e1/e2*hAir + e1/eBar*hBar);
	else
		E = gammaMat ./ betaMat .* (u2/u1*hAir + uBar/u1*hBar);
	end
end

% all parameters are matrices
% evaluates the following:
%		c1*e2*(kxn*cc*i + k*ss)
%		- c1*e1*(kxn*cc*i - k*ss)
%		+ c2*e2*(kxn*ss*i - k*cc)
%		- c2*e1*(-kxn*ss*i - k*cc)
function [val] = eval(kxn, k, c1, c2, e1, e2, cc, ss)
	t1 = kxn.*cc.*1i;
	t2 = kxn.*ss.*1i;
	t3 = k.*ss;
	t4 = k.*cc;
	
	val = c1.*(e2.*(t1+t3) - e1.*(t1-t3)) + c2.*(e2.*(t2-t4) - e1.*(-t2-t4));
end