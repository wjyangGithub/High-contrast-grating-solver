%% function [powRefl, powTrans, fldRefl, fldTrans, superModeSingleTripPhase] = solveRT(H, E, period, lambda, e1, u1, beta, gamma, tg, N, M, isTM)
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
% This function calculates reflectivity and transmission for various orders.

% Input parameters (H and E are [N x M] matrices, 
%                   beta is a [1 x M] vectors, gamma is a [1 x N] vector; 
%                   all others are scalars)
%	H			H field mode overlap integral				    [1]
%	E			E field mode overlap integral			        [1]
%   beta        z-direction mode propagation constants in the waveguide arrays   [1/m]
%   gamma       z-direction mode propagation constants outside the HCG           [1/m]
%	tg			thickness of grating							[m]
%	lambda		wavelength of light	in vacuum					[m]
%	e1			relative permittivity above and below grating	[1]
%	u1			relative permeability above and below grating	[1]
%	period		period of HCG									[m]
%	N			number of orders outside grating (total: N)	    [1]
%	M			number of orders inside grating	(total: M)		[1]
%	isTM		1 for TM polarization, 0 for TE					[boolean]

% Output parameters (powRefl, powTrans, fldRefl, and fldTrans are [1 X (2N+1)] vectors, 
%                   superModeSingleTripPhase is a [1 x M] matrix)
%   powRefl     power reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1]
%   powTrans    power transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1] 
%   fldRefl     field reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1] 
%   fldTrans    field transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1]
%   superModeSingleTripPhase       the single trip phase of the supermode inside the HCG  [rad]  

function [powRefl, powTrans, fldRefl, fldTrans, superModeSingleTripPhase] = solveRT(H, E, period, lambda, e1, u1, beta, gamma, tg, N, M, isTM)
	phi = diag(exp(-1i*beta*tg));	% phi is diagonal
	
	if (N == M)
		rho = (H+E) \ (H-E);
	else
		rho = pinv(H+E) * (H-E);
	end
	
	% precalculate
	prp = phi*rho*phi;
	Hprp = H*prp;
	Eprp = E*prp;
	
	if (N == M)
		R = (E+Eprp - H+Hprp) / (H-Hprp + E+Eprp);
		T = (2*E*phi + 2*E*rho*phi) / (H-Hprp + E+Eprp);
	else
		R = (E+Eprp - H+Hprp) * pinv(H-Hprp + E+Eprp);
		T = (2*E*phi + 2*E*rho*phi) * pinv(H-Hprp + E+Eprp);
	end
	
	if (isTM)
        gammaMult = gamma ./ gamma(1);
	else
		gammaMult = gamma(1) ./ gamma;
	end
		
	fldRefl = conj(transpose(R(:,1))) .* sqrt(gammaMult);	
	fldTrans = conj(transpose(T(:,1))).* sqrt(gammaMult);    
    fldRefl(2:end)=fldRefl(2:end)/2;
    fldTrans(2:end)=fldTrans(2:end)/2;
	
	% artificially set evanescent wave to 0
	thresh = floor(period * sqrt(e1*u1) / lambda) + 1;
	fldRefl(thresh+1:end) = 0;
	fldTrans(thresh+1:end) = 0;
	
	% mirror to make symmetric
	fldRefl = [fldRefl(end:-1:2), fldRefl];
	fldTrans = [fldTrans(end:-1:2), fldTrans];
	
    powRefl = abs(fldRefl).^2;	
	powTrans = abs(fldTrans).^2;
    
    % The superModeSingleTripPhase describes the single trip phase of the supermode inside the HCG.
    superModeSingleTripPhase=transpose(angle(eig(rho*phi)));
end