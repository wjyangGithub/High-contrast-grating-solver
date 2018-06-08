%% function [ka, ks, beta] = solveKaKsBeta_OA(a, s, lambda, eBar, uBar, e1, u1, e2, u2, theta, M, isTM)	
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
%	e1			relative permittivity above and below grating	[1]
%	u1			relative permeability above and below grating	[1]
%	e2			relative permittivity between grating bars		[1]
%	u2			relative permeability between grating bars		[1]
%	theta		incident light angle							[rad]
%	M			number of solutions to solve for				[1]
%	isTM		1 for TM polarization, 0 for TE					[boolean]

% Output parameters ([1 x M] vectors):
%	ka			transverse wave vectors in the low index air gaps        [1/m]
%	ks			transverse wave vectors in the high index grating bars	 [1/m]
%   beta        z-direction mode propagation constants in the waveguide arrays       [1/m]

function [ka, ks, beta] = solveKaKsBeta_OA(a, s, lambda, eBar, uBar, e1, u1, e2, u2, theta, M, isTM)	
	P = cos(2*pi*sqrt(e1*u1)*(a+s)/lambda*sin(theta));	% constant
	T = 2*pi/lambda * sqrt(eBar*uBar - e2*u2);		    % threshold = sqrt(U)
	U = T^2;
	
	if (isTM)
		N = e2 / eBar;
	else
		N = 1;
	end
	
	sols = zeros(1, M);
	
	n1Max = floor(s*T/pi - 1/2);
	n1 = 0:1:n1Max;
	
	% also search last extremum to threshold
	% in case no solutions, search from 1e-3 (to exclude 0) to threshold
	bounds1 = [1e-3, rootBound1(n1, s), T-1e-3];
	
	found = 0;		% keep track of number of solutions found
	for i=1:length(bounds1)-1
		curr = searchInterval(bounds1(i), bounds1(i+1), a, s, P, N, U, false);
		sols(found+1:1:found+length(curr)) = curr;
		found = found + length(curr);
		
		if (found >= M)		% if done, stop
			ks = sols(1:M);
			ka = sqrt(ks.^2 - U);
			return;
		end
	end
	
	% bounds for roots past threshold
	nMin2 = ceil(s*T/pi);
	
	temp = rootBound2(nMin2, a, s, U);	% since ks is not strictly increasing when a ~= s, ensure past threshold
	if (temp <= T + 2)
		nMin2 = nMin2 + 1;
	end
	
	nMax2 = nMin2 + (M-found);
	n2 = nMin2:1:nMax2;
	bounds2 = rootBound2(n2, a, s, U);
	
	% search threshold to first extrema
	errMin = err(T+1e-3, a, s, P, N, U, true);
	errMax = err(bounds2(1), a, s, P, N, U, true);
	test = findRoot(T+1e-3, bounds2(1), a, s, P, N, U, true, errMin, errMax);
	if (test ~= -1)
		found = found + 1;
		sols(found) = test;
	end
	
	% stop if done
	if (found == M)
		ks = sols;
		ka = sqrt(ks.^2 - U);
		return;
	end
	
	% find remaining necessary solutions
	errMin = errMax;
	for i=1:M-found
		errMax = err(bounds2(i+1), a, s, P, N, U, true);
		sols(i+found) = findRoot(bounds2(i), bounds2(i+1), a, s, P, N, U, true, errMin, errMax);
		errMin = errMax;
	end
	
	ks = sols;
	ka = sqrt(ks.^2 - U);
    
    beta = conj(sqrt(eBar*uBar*(2*pi/lambda).^2-ks.^2));

end

% bounds for roots for ks < T
function [ks] = rootBound1(n, s)
	ks = pi/(2*s) + n.*pi/s;
end

% bounds for roots for ks > T
function [ks] = rootBound2(n, a, s, U)
	if (a == s)
		ks = (n.^2*pi^2 + a^2*U) ./ (2*n.*pi*s);
	else
		ks = (n.*pi*s - a*sqrt(n.^2*pi^2 + U*(a^2 - s^2))) / (s^2-a^2);
	end
end

% bisection method for roots
% don't need to recalculate err at one endpoint each iteration
function [root] = findRoot(min, max, a, s, P, N, U, realKa, errMin, errMax)
	%fprintf('%f\t%f\n', min, max);
	if ((errMin>0) == (errMax>0))		% no sign change
		root = -1;
		%fprintf('no solution\n\n');
		return;
	end
	
	mid = (min+max) / 2;
	if (max - min < 1e2)		% set precision here
		root = mid;
		%fprintf('%f\n\n', mid);
		return;
	end
	
	errMid = err(mid, a, s, P, N, U, realKa);
	if ((errMin>0) == (errMid>0))
		root = findRoot(mid, max, a, s, P, N, U, realKa, errMid, errMax);
	else
		root = findRoot(min, mid, a, s, P, N, U, realKa, errMin, errMid);
	end
end

% determine whether 0, 1, or 2 sols in this interval
function [sols] = searchInterval(min, max, a, s, P, N, U, realKa)
	%fprintf('searching interval %f to %f\n', min, max);
	errMin = err(min, a, s, P, N, U, realKa);
	errMax = err(max, a, s, P, N, U, realKa);
	
	% no sign change; see if 0 or 2 solutions
	if ((errMin>0) == (errMax>0))
		derivMin = err(min+1, a, s, P, N, U, realKa) - errMin;
		derivMax = err(max+1, a, s, P, N, U, realKa) - errMax;
		%fprintf('no sign change: dMin = %f\t\tdMax = %f\n', derivMin, derivMax);
		%fprintf('\t\teMin = %f\t\teMax = %f\n', errMin, errMax);
		
		% 2 solutions
		if (((derivMin>0) ~= (derivMax>0)) && ((derivMin>0) ~= (errMin>0)))
			% find a value with the opposite sign between min and max
			mid = findMidVal(min, max, derivMin>0, a, s, P, N, U, realKa, derivMin);
			errMid = err(mid, a, s, P, N, U, realKa);
			%fprintf('2 solutions: midVal = %f\n', mid);
			
			sols = zeros(1,2);
			sols(1) = findRoot(min, mid, a, s, P, N, U, realKa, errMin, errMid);
			sols(2) = findRoot(mid, max, a, s, P, N, U, realKa, errMid, errMax);
		
		% no solutions; return empty
		else
			%fprintf('0 solutions\n');
			sols = zeros(0,1);
		end
		
	% there is a sign change; use bisection	
	else
		%fprintf('1 solution\n');
		sols = findRoot(min, max, a, s, P, N, U, realKa, errMin, errMax);
	end
end

% derivative of error function
function [deriv] = errorDeriv(ks, a, s, P, N, U, realKa)
	step = 1;
	d1 = err(ks, a, s, P, N, U, realKa);
	d2 = err(ks+step, a, s, P, N, U, realKa);
	
	deriv = (d2-d1) / step;
end

% bisection to find value with a certain sign
% searchPos: looking for positive value or negative value
function [midVal] = findMidVal(min, max, searchPos, a, s, P, N, U, realKa, derivMin)
	mid = (min+max) / 2;
	errMid = err(mid, a, s, P, N, U, realKa);
	if ((errMid > 0) == searchPos)
		midVal = mid;
		return;
	elseif (max - min < 1)
		%fprintf('Error: no opposite sign\n');
		return;
	end
	
	derivMid = errorDeriv(mid, a, s, P, N, U, realKa);
	if ((derivMin>0) == (derivMid>0))
		midVal = findMidVal(mid, max, searchPos, a, s, P, N, U, realKa, derivMid);
	else
		midVal = findMidVal(min, mid, searchPos, a, s, P, N, U, realKa, derivMin);
	end
end

% The error function - value is always real
function errVal = err(ks, a, s, P, N, U, realKa)	
	% use formula for real ka
	if (realKa)
		ka = sqrt(ks^2 - U);
		errVal = 2*N*ka*ks * (P - cos(ka*a)*cos(ks*s)) + (ka^2 + N^2*ks^2) * sin(ka*a)*sin(ks*s);
		
	% use formula for imag ka
	else
		% precalculating a little
		G = sqrt(U - ks^2);					% actual ka / i
		B = 2*N*G*ks*cosh(G*a);
		C = (G^2-N^2*ks^2) * sinh(G*a);
		
		% normalize to change derivative properties, but zeros are still zeros
		errVal = (2*N*G*ks*P - B*cos(ks*s) - C*sin(ks*s)) / sqrt(B^2+C^2);
	end
end