%% function [powRefl, powTrans, fldRefl, fldTrans, superModeSingleTripPhase] = solveRT_OA(H, E, beta, gamma, tg, lambda, e1, u1, period, theta, N, M, isTM)
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


% This function calculates reflectivity and transmission for various orders.

% Input parameters (H and E are [(2N+1) x M] matrices, 
%                   beta is a [1 x M] vectors, gamma is a [1 x (2N+1)] vector; 
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
%	theta		incident light angle							[rad]
%	N			number of orders outside grating (total: 2N+1)	[1]
%	M			number of orders inside grating	(total: M)		[1]
%	isTM		1 for TM polarization, 0 for TE					[boolean]

% Output parameters (powRefl, powTrans, fldRefl, and fldTrans are [1 X (2N+1)] vectors, 
%                   superModeSingleTripPhase is a [1 x M] matrix)
%   powRefl     power reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1]
%   powTrans    power transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1] 
%   fldRefl     field reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1] 
%   fldTrans    field transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1]
%   superModeSingleTripPhase       the single trip phase of the supermode inside the HCG  [rad]  

function [powRefl, powTrans, fldRefl, fldTrans, superModeSingleTripPhase] = solveRT_OA(H, E, beta, gamma, tg, lambda, e1, u1, period, theta, N, M, isTM)
	phi = diag(exp(-1i*beta*tg), 0);	% phi is a diagonal matrix
	
	if (2*N+1 == M)
		rho = (H+E) \ (H-E);
	else
		rho = pinv(H+E) * (H-E);
	end
	
	% precalculate
	prp = phi*rho*phi;
	Hprp = H*prp;
	Eprp = E*prp;
	
	if (2*N+1 == M)
		t1 = H-Hprp + E+Eprp;
		R = (E+Eprp - H+Hprp) / t1;
		T = (2*E*phi + 2*E*rho*phi) / t1;
	else
		invMat = pinv(H-Hprp + E+Eprp);
		R = (E+Eprp-H+Hprp) * invMat;
		T = (2*E*phi + 2*E*rho*phi) * invMat;
    end
		
    fldRefl = conj(transpose(R(:,N+1)));
    fldTrans = conj(transpose(T(:,N+1)));
        
	% artificially set evanescent orders to 0
	n1 = sqrt(e1*u1);
	minN = floor(n1*period/lambda*(-1-sin(theta))) + N + 1;
	maxN = ceil(n1*period/lambda*(1-sin(theta))) + N + 1;
	
	fldRefl([1:1:minN, maxN:1:2*N+1]) = 0;
	fldTrans([1:1:minN, maxN:1:2*N+1]) = 0;
	
	% correct orders other than 0th order; 0's are still 0
	if (isTM)
        gammaMult = gamma ./ gamma(N+1);
	else
        gammaMult = gamma(N+1) ./ gamma;
	end
	
	fldRefl = fldRefl .* sqrt(gammaMult);
	fldTrans = fldTrans .* sqrt(gammaMult);
    
    powRefl = abs(fldRefl).^2;	
	powTrans = abs(fldTrans).^2;
    
    % The superModeSingleTripPhase describes the single trip phase of the supermode inside the HCG.
    superModeSingleTripPhase=transpose(angle(eig(rho*phi)));

end