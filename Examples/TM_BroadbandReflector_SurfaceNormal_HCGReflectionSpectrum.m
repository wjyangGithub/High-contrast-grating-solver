%% TM_BroadbandReflector_SurfaceNormal_HCGReflectionSpectrum.m
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
% This program calculates the HCG reflection and transmission spectrum with a plane wave
% incident at surface normal incident angle.

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

clc;
clear all;
close all;

currentFolder = pwd;
cd ../;
cd('./SurfaceNormalIncidence');

%% HCG parameters
% User defined parameters
period=0.77e-6;             % [m] HCG period
DC=0.76;                    % [1] HCG duty cycle
tg=0.455e-6;                % [m] HCG thickness
s=DC*period;                % [m] grating bar width
a=(1-DC)*period;            % [m] air gap width
nbar=3.48;                  % [1] refractive index of grating bar
eBar=nbar^2;                % [1] relative permittivity of grating bar
uBar=1;                     % [1] relative permeability of grating bar
e1=1;                       % [1] relative permittivity above and below grating
u1=1;                       % [1] relative permeability above and below grating
e2=1;                       % [1] relative permittivity of the grating gap
u2=1;                       % [1] relative permeability of the grating gap
ord=10;                     % [1] orders
isTM=1;                     % [boolean] 1 for TM, 0 for TE

% set the orders for N and M. M, N is set to be the same, so that the matrix
% involved in the calculation is square.
N = ord;
M = ord;

% set the wavelength to be scanned
lambda=(1.0:0.001:3.0)*1e-6;        % [m] incident wave wavelength
lambda=lambda+(1e-6)*1e-6;          % [m] incident wave wavelength, so as to avoid gamma=0

%% Calculation of the reflection spectrum
powRefl0th=zeros(1,length(lambda));    % [1] 0th order power reflection coefficient
powTrans0th=zeros(1,length(lambda));   % [1] 0th order power transmission coefficient
fldRefl0th=zeros(1,length(lambda));    % [1] 0th order field reflection coefficient
fldTrans0th=zeros(1,length(lambda));   % [1] 0th order field transmission coefficient

counter=1;
h = waitbar(0,'Calculating...');
tic;
for idx=1:length(lambda)
    [ka, ks, beta] = solveKaKsBeta(a, s, eBar, uBar, e2, u2, lambda(idx), ord+1, isTM);
	[H, E, ~, gamma] = solveHE(ka, ks, a, s, eBar, uBar, e1, u1, e2, u2, a+s, lambda(idx), ord+1, ord+1, isTM);
	[powRefl, powTrans, fldRefl, fldTrans,~] = solveRT(H, E, a+s, lambda(idx), e1, u1, beta, gamma, tg, ord+1, ord+1, isTM);
    powRefl0th(idx)=powRefl((1+end)/2);         % extract the 0th order power reflection coefficient
    powTrans0th(idx)=powTrans((1+end)/2);       % extract the 0th order power transmission coefficient
    fldRefl0th(idx)=fldRefl((1+end)/2);         % extract the 0th order field reflection coefficient
    fldTrans0th(idx)=fldTrans((1+end)/2);       % extract the 0th order field transmision coefficient

%%%% Alternatively, use "solveHCG.m" in the following code    
%{
    [fldRefl, fldTrans] = solveHCG(a, s, eBar, uBar, e1, u1, e2, u2, lambda(idx), tg, ord, isTM);
    fldRefl0th(idx)=fldRefl((1+end)/2);         % extract the 0th order field reflection coefficient
    fldTrans0th(idx)=fldTrans((1+end)/2);       % extract the 0th order field transmision coefficient
    powRefl0th(idx)=(abs(fldRefl0th(idx)))^2;   % the 0th order power reflection coefficient
    powTrans0th(idx)=(abs(fldTrans0th(idx)))^2; % the 0th order power transmission coefficient
%}

    if (idx>length(lambda)/10*counter)
        waitbar(idx/length(lambda),h);
        counter=counter+1;
    end
end
close(h);
tElapsed=toc;
fprintf(['Total Calculation Time = ' num2str(tElapsed) ' seconds \n']);
%   Note N=ord
%   powRefl     power reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1]
%   powTrans    power transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1] 
%   fldRefl     field reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1] 
%   fldTrans    field transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1]

%% Plot the figure
figure;  % Power reflectivity
plot(lambda*1e6,powRefl0th,'r');
xlabel('Wavelength (\mum)');
ylabel('Power Reflectivity');

figure;  % Power transmission coefficient
plot(lambda*1e6,powTrans0th,'r');
xlabel('Wavelength (\mum)');
ylabel('Power Transmission Coefficient');

figure;  % Reflection phase
plot(lambda*1e6,angle(fldRefl0th),'r');
xlabel('Wavelength (\mum)');
ylabel('Reflection Phase (rad)');

figure;  % Transmission phase
plot(lambda*1e6,angle(fldTrans0th),'r');
xlabel('Wavelength (\mum)');
ylabel('Transmission Phase (rad)');

%% Return to original folder
cd (currentFolder);
