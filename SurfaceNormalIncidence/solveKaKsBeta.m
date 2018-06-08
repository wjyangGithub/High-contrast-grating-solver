%% function [ka, ks, beta] = solveKaKsBeta(a, s, eBar, uBar, e2, u2, lambda, M, isTM)
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
% This function solves the eigen modes in waveguide arrays. It returns the
% propagation constant beta and the transverse k-vectors in the high index
% grating bars and the low index air gap.

% This function returns the first M solutions for ka, ks and beta, given that ka, ks ~= 0
% This function solves the transcendental equation using a semi-analytic algorithm. 
% Bounds for the solutions are determined analytically, while the bisection method is employed
% to refine each solution to an arbitrary precision.

% Input parameters (all parameters are scalars):
%	a			air gap	width									[m]
%	s			grating bar width								[m]
%	lambda		wavelength of light	in vacuum					[m]
%	eBar		relative permittivity of grating bar			[1]
%	uBar		relative permeability of grating bar			[1]
%	e2			relative permittivity between grating bars		[1]
%	u2			relative permeability between grating bars		[1]
%	M			number of solutions to solve for				[1]
%	isTM		1 for TM polarization, 0 for TE					[boolean]

% Output parameters ([1 x M] vectors):
%	ka			transverse wave vectors in the low index air gaps        [1/m]
%	ks			transverse wave vectors in the high index grating bars	 [1/m]
%   beta        z-direction mode propagation constants in the waveguide arrays       [1/m]

function [ka, ks, beta] = solveKaKsBeta(a, s, eBar, uBar, e2, u2, lambda, M, isTM)
	% define widely used constants
	U = (2*pi/lambda)^2 * (eBar*uBar - e2*u2);
	if (isTM)
		N = e2 / eBar;
	else
		N = 1;
	end
	
	% find M+1 vertical asymptotes to guarantee at least M solutions
	asym = asymptotes(s, a, U, M);
	
	% initialize array
	ks = zeros(1,M);
	
	% check whether there is a solution between 0 and first asymptote
	% noninclusive; if there is, only consider the first M asymptotes
	iOff = 0;
	iEnd = M;
	errMin = error(1e-5, s, a, N, U);
	errMax = error(asym(1)-1e-5, s, a, N, U);
	sol0 = findSolution(s, a, N, U, 1e-5, asym(1)-1e-5, errMin, errMax);
	if (sol0 ~= -1)
		ks(1) = sol0;
		iOff = 1;
		iEnd = M - 1;
	end
	
	% find solutions between adjacent asymptotes
	for i=1:1:iEnd
		errMin = error(asym(i)+1e-5, s, a, N, U);
		errMax = error(asym(1+1)-1e-5, s, a, N, U);
		ks(i+iOff) = findSolution(s, a, N, U, asym(i)+1e-5, asym(i+1)-1e-5, errMin, errMax);
	end
	
	% find ka given ks
	ka = toKa(ks, U);
    
    % find beta given ks
    beta = conj( sqrt((2*pi*sqrt(eBar*uBar)/lambda)^2 - ks.^2) );		% (1 x M vector)

end

% Recursively find the root to error function between
% min and max to a precision +/- prec
function sol = findSolution(s, a, N, U, min, max, errMin, errMax)
	if ((errMin>0) == (errMax>0))	% no solution if errMin has same sign as errMax
		sol = -1;
		return;
	end
	
	mid = (min + max) / 2;
	if ((max - min) < 1)				% precision up to 1e-8
		sol = mid;
		return;
	end	
	
	errMid = error(mid, s, a, N, U);
	if ((errMin>0) == (errMid>0))				% min and mid have same sign; search (mid, max)
		sol = findSolution(s, a, N, U, mid, max, errMid, errMax);
	else										% max and mid have same sign; search (min, mid)
		sol = findSolution(s, a, N, U, min, mid, errMin, errMid);
	end
end

% Find the asymptotes of the error function within range [min, max]
function [asym] = asymptotes(s, a, U, M)
	% loop one by one (faster than vectorizing then sorting)
	n1 = -1;
	n2 = -1;
	
	asym = zeros(1, M+1);
	
	k1 = ks1(0, s);
	k2 = ks2(0, a, U);
	for m = 1:M+1
		if (abs(k1 - k2) < 1e-5)		% will happen if TE and a/s is a ratio of positive odd integers <= M
			asym(m) = k1;
			n1 = n1 + 1;
			k1 = ks1(n1+1, s);
			n2 = n2 + 1;
			k2 = ks2(n2+1, a, U);
		elseif (k1 < k2)
			asym(m) = k1;
			n1 = n1 + 1;
			k1 = ks1(n1+1, s);
		else
			asym(m) = k2;
			n2 = n2 + 1;
			k2 = ks2(n2+1, a, U);
		end
	end
end

% Asymptotes for the first tangent of the error function
% n;			% vector of integers
function [ks] = ks1(n, s)
	ks = 2*pi*n/s + pi/s;
end

% Asymptotes for the second tangent of the error function
% n;			% vector of integers
function [ks] = ks2(n, a, U)
	ks = sqrt((2*pi*n/a + pi/a).^2 + U);
end

% One variable function of ks that we are trying to find the roots of
function err = error(ks, s, a, N, U)
	ka = toKa(ks, U);
	err = N * ks*tan(ks * s/2) + ka*tan(ka * a/2);
end

% ka = sqrt( ks^2 - (2*pi/lambda)^2 * (nBar^2-1) )
function [ka] = toKa(ks, U)
	ka = sqrt(ks.^2 - U);
end