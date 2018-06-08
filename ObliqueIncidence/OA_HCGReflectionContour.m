%% OA_HCGReflectionContour.m
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


% This program calculates the HCG reflectivity contour against HCG thickness 
% and incident light wavelength, with a plane wave incident at an incident angle.

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
%	theta		incident light angle							[rad]
%	ord			number of orders								[1]
%	isTM		1 for TM, 0 for TE								[boolean]

clc;
clear all;
close all;

%% HCG parameters
% User defined parameters
period=1e-6;                % [m] HCG period
DC=0.6;                     % [1] HCG duty cycle
s=DC*period;                % [m] grating bar width
a=(1-DC)*period;            % [m] air gap width
nbar=3.48;                  % [1] refractive index of semiconductor
eBar=nbar^2;                % [1] relative permittivity of grating bar
uBar=1;                     % [1] relative permeability of grating bar
e1=1;                       % [1] relative permittivity above and below grating
u1=1;                       % [1] relative permeability above and below grating
e2=1;                       % [1] relative permittivity of the grating gap
u2=1;                       % [1] relative permeability of the grating gap
theta=50/180*pi;            % [rad] incident angle, with respect to the normal
ord=8;                      % [1] orders
isTM=1;                     % [boolean] 1 for TM, 0 for TE

% set the orders for N and M. M is set to be 2N+1, so that the matrix
% involved in the calculation is square.
N = ord;
M = 2*N + 1;

% set the wavelength and thickness to be scanned
lambda=(1:0.02:4)*period;         % [m] incident wave wavelength
lambda=lambda+(1e-6)*period;      % [m] incident wave wavelength, so as to avoid gamma=0
tg=(0:0.02:3)*period;             % [m] HCG thickness

%% Calculation of reflection contour
powRefl0th=zeros(length(tg),length(lambda));
powTrans0th=zeros(length(tg),length(lambda));

counter=1;
h = waitbar(0,'Calculating...');
tic;
for idx2=1:length(lambda)
    [ka, ks, beta] = solveKaKsBeta_OA(a, s, lambda(idx2), eBar, uBar, e1, u1, e2, u2, theta, M, isTM);
    [H, E, ~, gamma] = solveHE_OA(ka, ks, a, s, a+s, eBar, uBar, e1, u1, e2, u2, lambda(idx2), theta, N, M, isTM);
    for idx1=1:length(tg)
        [powRefl, powTrans, fldRefl, fldTrans, ~] = solveRT_OA(H, E, beta, gamma, tg(idx1), lambda(idx2), e1, u1, a+s, theta, N, M, isTM);
        powRefl0th(idx1,idx2)=powRefl((1+end)/2);        % extract the 0th order power reflection coefficient
        powTrans0th(idx1,idx2)=powTrans((1+end)/2);      % extract the 0th order power transmission coefficient
    end
    
    if (idx2>length(lambda)/10*counter)
        waitbar(idx2/length(lambda),h);
        counter=counter+1;
    end
end
close(h);
tElapsed=toc;
fprintf(['Total Calculation Time = ' num2str(tElapsed) ' seconds \n']);
fprintf(['Plotting figures... \n']);
%   Note N = ord;
%   powRefl     power reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1]
%   powTrans    power transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1] 

%% Plot reflectivity contour
figure;
v=0:0.02:1;
contourf(lambda./period, tg./period, powRefl0th, v, 'linestyle','none');
xlabel('\lambda/\Lambda');
ylabel('tg/\Lambda');
h=colorbar;
set(get(h,'ylabel'),'String','Reflectivity','Rotation',270,'fontSize',26,'fontWeight','bold');
h1=get(get(h,'ylabel'),'position');
h1(1)=h1(1)+5;
set(get(h,'ylabel'),'position',h1);
caxis([0 1]);
axis([1 4 0 3]);