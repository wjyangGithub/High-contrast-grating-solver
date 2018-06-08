%% function [fldRefl, fldTrans] = solveHCG(a, s, eBar, uBar, e1, u1, e2, u2, lambda, tg, ord, isTM)
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
% This program calculates the HCG field reflection and transmission coefficient
% with a plane wave incident at surface normal incident angle.

% Input parameters (all parameters are scalars):
%	a			air gap	width									[m]
%	s			grating bar width								[m]
%	lambda		incident light wavelength in vacuum				[m]
%	tg			grating thickness								[m]
%	eBar		relative permittivity of grating bar			[1]
%	uBar		relative permeability of grating bar			[1]
%	e1			relative permittivity above and below grating	[1]
%	u1			relative permeability above and below grating	[1]
%	e2			relative permittivity between grating bars		[1]
%	u2			relative permeability between grating bars		[1]
%	ord			number of orders								[1]
%	isTM		1 for TM, 0 for TE								[boolean]

% Output parameters (fldRefl, and fldTrans are [1 X (2N+1)] vectors, N=ord+1
%   fldRefl     field reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1] 
%   fldTrans    field transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1]

% Note	
%	6 orders is enough for subwavelength gratings.
%	20 orders is needed for non-subwavelength gratings.
%   lambda ~= (a+s)	->	shift lambda very slightly e.g. 1e-12


function [fldRefl, fldTrans] = solveHCG(a, s, eBar, uBar, e1, u1, e2, u2, lambda, tg, ord, isTM)
	warning('off', 'MATLAB:nearlySingularMatrix');
	
	[ka, ks, beta] = solveKaKsBeta(a, s, eBar, uBar, e2, u2, lambda, ord+1, isTM);
	[H, E, ~, gamma] = solveHE(ka, ks, a, s, eBar, uBar, e1, u1, e2, u2, a+s, lambda, ord+1, ord+1, isTM);
	[~, ~, fldRefl, fldTrans, ~] = solveRT(H, E, a+s, lambda, e1, u1, beta, gamma, tg, ord+1, ord+1, isTM);
	
	warning('on', 'MATLAB:nearlySingularMatrix');
end