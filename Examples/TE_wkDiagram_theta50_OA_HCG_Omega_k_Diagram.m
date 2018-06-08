%% TE_wkDiagram_theta50_OA_HCG_Omega_k_Diagram.m
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


% This program calculates the omega-k diagram of the waveguide array
% (region II of HCG) with a plane wave incident at an incident angle.
% This program calculates both the even order modes and odd order modes,
% though only even order modes are excited at surface normal incidence.

% Input parameters (all parameters are scalars):
%	a			air gap	width								    [m]
%	s			grating bar width								[m]
%	lambda		incident light wavelength in vaccum				[m]
%	eBar		relative permittivity of grating bar			[1]
%	uBar		relative permeability of grating bar			[1]
%	e1			relative permittivity above and below grating	[1]
%	u1			relative permeability above and below grating	[1]
%	e2			relative permittivity between grating bars		[1]
%	u2			relative permeability between grating bars		[1]
%	theta		incident light angle							[rad]
%   M           number of modes to search                       [1]
%	isTM		1 for TM, 0 for TE								[boolean]

clc;
clear all;
close all;

currentFolder = pwd;
cd ../;
cd('./ObliqueIncidence');

%% General parameters
c=299792458;                % [m/s] speed of light in vaccum

%% HCG parameters
period=1e-6;                % [m] HCG period
DC=0.5;                     % [1] HCG duty cycle
s=DC*period;                % [m] grating bar width
a=(1-DC)*period;            % [m] air gap width
nbar=3.48;                  % [1] refractive index of grating bar
n2=1;                       % [1] refractive index of air gap
eBar=nbar^2;                % [1] relative permittivity of grating bar
uBar=1;                     % [1] relative permeability of grating bar
e1=1;                       % [1] relative permittivity above and below grating
u1=1;                       % [1] relative permeability above and below grating
e2=n2^2;                    % [1] relative permittivity of the grating gap
u2=1;                       % [1] relative permeability of the grating gap
theta=50/180*pi;             % [rad] incident angle, with respect to the normal
M=10;                       % [1] search for M modes
isTM=0;                     % [boolean] 1 for TM, 0 for TE

omegaSet=0+pi*c/period/50:pi*c/period/401:4*pi*c/period;        % [rad/s] incident wave angular frequency, in rad/s
lambdaSet=2*pi*c./omegaSet;                                     % [m] incident wave wavelength in vaccum
%  lambdaSet=(1:0.005:4)*period;                                % [m] incident wave wavelength in vaccum
%  lambdaSet=lambdaSet+(1e-6)*period;                           % [m] incident wave wavelength in vaccum
%  omegaSet=2*pi*c./lambdaSet;                                  % [rad/s] incident wave frequency, in rad/s

%% Calculation of omege-k diagram
betaSet=zeros(length(lambdaSet),M);
ksSet=zeros(length(lambdaSet),M);
kaSet=zeros(length(lambdaSet),M);

counter=1;
h = waitbar(0,'Calculating...');
tic;
for idx=1:length(lambdaSet)
    [ka, ks, beta] = solveKaKsBeta_OA(a, s, lambdaSet(idx), eBar, uBar, e1, u1, e2, u2, theta, M, isTM);
    kaSet(idx,:)=ka;
    ksSet(idx,:)=ks;
    betaSet(idx,:)=beta;
    if (idx>length(lambdaSet)/10*counter)
        waitbar(idx/length(lambdaSet),h);
        counter=counter+1;
    end
end
close(h);
tElapsed=toc;
fprintf(['Total Calculation Time = ' num2str(tElapsed) ' seconds \n']);

%% Plot omega-k diagram
figure;
hold on;
plot(omegaSet/(2*pi*c/period),omegaSet*nbar/((2*pi/period)*c),'k--', 'lineWidth', 6);
plot(omegaSet/(2*pi*c/period),omegaSet*n2/((2*pi/period)*c),'k--', 'lineWidth', 6);
plot(omegaSet/(2*pi*c/period), real(betaSet/(2*pi/period)), 'lineWidth', 6);
xlabel('\omega [units of 2\pic/\Lambda]','fontsize',30);
ylabel('\beta [units of 2\pi/\Lambda]','fontsize',30);
set(gca,'LineWidth',4, 'Fontsize',30, 'Ytick', 0:2:10);
axis([0 2 0.02 6]);
box; 

%% Return to original folder
cd (currentFolder);