%% HCGFieldProfile.m
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
% This program calculates the electromagnetic field profile in the HCG
% with a plane incident at surface normal incident angle.

% Input parameters (all parameters are scalars):
%   period      HCG period                                      [m]
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
%   x;          x vector to mesh the HCG in the x direction     [m]
%   z;          z vector to mesh the HCG in the z direction     [m]

clc;
clear all;
close all;

%% HCG parameters
% User defined parameters
period=1.12e-6;             % [m] HCG period
DC=0.4;                     % [1] HCG duty cycle
tg=0.52e-6;                 % [m] HCG thickness
s=DC*period;                % [m] grating bar width
a=(1-DC)*period;            % [m] air gap width
nbar=3.48;                  % [1] refractive index of semiconductor
eBar=nbar^2;                % [1] relative permittivity of grating bar
uBar=1;                     % [1] relative permeability of grating bar
e1=1;                       % [1] relative permittivity above and below grating
u1=1;                       % [1] relative permeability above and below grating
e2=1;                       % [1] relative permittivity of the grating gap
u2=1;                       % [1] relative permeability of the grating gap
ord=10;                     % [1] orders
isTM=0;                     % [boolean] 1 for TM, 0 for TE

lambda=1.55e-6;                          % [m] incident light wavelength in vacuum
plotLength=1;                             % [1] indicator of the length of the plot in z direction    
x=(0:1e-3:1)*period;                      % [m] mesh grid for x
z=(-2*tg-plotLength*tg):(period*1e-3):(tg+plotLength*tg); % [m] mesh grid for z

%% Calculation of field profile
tic;
[~, ~, H_all, E_all, incidentWaveAmp ] = solveFieldProfile( a, s, eBar, uBar, e1, u1, e2, u2, lambda, tg, ord, isTM, x, z );
tElapsed=toc;
fprintf(['Total Calculation Time = ' num2str(tElapsed) ' seconds \n']);
fprintf('Plotting Figures...\n');

%% Plot the field profile: E field  (TM: Ex; TE: Ey)
Z_E_intensity=[abs(E_all).^2/(incidentWaveAmp(2))^2 ...
    abs(E_all).^2/(incidentWaveAmp(2))^2 ...
    abs(E_all).^2/(incidentWaveAmp(2))^2];
Z_E_field=[E_all./incidentWaveAmp(2) ...
    E_all./incidentWaveAmp(2) ...
    E_all./incidentWaveAmp(2)];
Z_E_field=real(Z_E_field);

figure();    % plot E field intensity
v=0:max(max(Z_E_intensity))/50:max(max(Z_E_intensity));
contourf([x./period-1 x./period x./period+1],fliplr(z./period), Z_E_intensity, v, 'linestyle','none' );
xlabel('x/\Lambda');
ylabel('z/\Lambda');
colorbar;
rectangle('Position',[a/period,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period-1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period+1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
daspect([1 1 1]);

figure();     % plot E field real part
v=min(min(Z_E_field)):(max(max(Z_E_field))-min(min(Z_E_field)))/50:max(max(Z_E_field));
contourf([x./period-1 x./period x./period+1],fliplr(z./period), Z_E_field, v, 'linestyle','none' );
xlabel('x/\Lambda');
ylabel('z/\Lambda');
colorbar;
rectangle('Position',[a/period,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period-1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period+1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
daspect([1 1 1]);


%% Plot the field profile: H field  (TM: Hy; TE: Hx)
Z_H_intensity=[abs(H_all).^2/(incidentWaveAmp(1))^2 ...
    abs(H_all).^2/(incidentWaveAmp(1))^2 ...
    abs(H_all).^2/(incidentWaveAmp(1))^2];
Z_H_field=[H_all./incidentWaveAmp(1) ...
    H_all./incidentWaveAmp(1) ...
    H_all./incidentWaveAmp(1)];
Z_H_field=real(Z_H_field);

figure();      % plot H field intensity
v=0:max(max(Z_H_intensity))/50:max(max(Z_H_intensity));
contourf([x./period-1 x./period x./period+1],fliplr(z./period), Z_H_intensity, v, 'linestyle','none' );
xlabel('x/\Lambda');
ylabel('z/\Lambda');
colorbar;
rectangle('Position',[a/period,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period-1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period+1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
daspect([1 1 1]);

figure();       % plot H field real part
v=min(min(Z_H_field)):(max(max(Z_H_field))-min(min(Z_H_field)))/50:max(max(Z_H_field));
contourf([x./period-1 x./period x./period+1],fliplr(z./period), Z_H_field, v, 'linestyle','none' );
xlabel('x/\Lambda');
ylabel('z/\Lambda');
colorbar;
rectangle('Position',[a/period,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period-1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
rectangle('Position',[a/period+1,-1*tg/period,s/period,tg/period], 'LineWidth',2,'LineStyle','--', 'EdgeColor',[1 1 1]);
daspect([1 1 1]);
