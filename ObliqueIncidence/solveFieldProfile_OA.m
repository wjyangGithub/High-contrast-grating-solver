%% function [ powRefl, powTrans, H_all, E_all, incidentWaveAmp ] = solveFieldProfile_OA( s, a, tg, eBar, uBar, e1, u1, e2, u2, lambda, theta, ord, isTM, x, z )
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


% This program calculates the electromagnetic field profile in the HCG
% with a plane wave incident at an incident angle.

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
%   x;          x vector to mesh the HCG in the x direction     [m]
%   z;          z vector to mesh the HCG in the z direction     [m]

% Output parameters (powRefl and powTrans are [1 X (2N+1)] vectors,
%                    H_all and E_all are matrixes with size defined by z and x: [length(z) X length(x)],
%                    incidentWaveAmp is 1 x 2 vectors)  
%   powRefl     power reflection coefficient for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)      [1]
%   powTrans    power transmission coefficent for various orders, formated in (-N, -(N-1), ...-2, -1, 0, 1, 2...N)     [1] 
%   H_all       H field profile (TM: Hy; TE: Hx)
%   E_all       E field profile (TM: Ex; TE: Ey)
%   incidentWaveAmp   amplitude of incident H field and E field 
%                     [incident H field, incident E field], (TM: Hy, Ex; TE: Hx, Ey) 

function [ powRefl, powTrans, H_all, E_all, incidentWaveAmp ] = solveFieldProfile_OA( s, a, tg, eBar, uBar, e1, u1, e2, u2, lambda, theta, ord, isTM, x, z )

period=a+s;             % [m]  HCG period    

% set the orders for N and M. M is set to be 2N+1, so that the matrix
% involved in the calculation is square.
N = ord;
M = 2*N + 1;

% set up parameters
e0=8.8542e-12;             % [F/m] vacuum permittivity
u0=4*pi*1e-7;              % [N*A^-2] vacuum permeability
eta0=sqrt(u0/e0);          % [ohm] vacuum wave impedance

n1 = sqrt(e1*u1);
k0 = 2*pi/lambda;            
kx0 = n1*k0 * sin(theta);
kz0 = n1*k0 * cos(theta);
n = -N:1:N;
kxn = kx0 + 2*pi*n/period;

%% Calculate the eigen modes in different regions
[ka, ks, beta] = solveKaKsBeta_OA(a, s, lambda, eBar, uBar, e1, u1, e2, u2, theta, M, isTM);

[H, E, ~, gamma] = solveHE_OA(ka, ks, a, s, a+s, eBar, uBar, e1, u1, e2, u2, lambda, theta, N, M, isTM);
index=find(x>a, 1, 'first');
deltaN0=zeros(1,2*N+1);
deltaN0((1+end)/2)=1;

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

%%%% TM polarization
if isTM
    hyin=zeros(M,length(x));
    exin=zeros(M,length(x));
    hyout=zeros(2*N+1,length(x));
    exout=zeros(2*N+1,length(x));
    	       
% matrix: hyin, exin, hyout, exout
    for m=1:1:M       
        hyin(m,1:index-1)=A(m)*cos(ka(m)*(x(1:index-1)-a/2))+B(m)*sin(ka(m)*(x(1:index-1)-a/2));
        hyin(m,index:end)=C(m)*cos(ks(m)*(x(index:end)-(a+period)/2))+D(m)*sin(ks(m)*(x(index:end)-(a+period)/2));        
        exin(m,1:index-1)=beta(m)/k0*eta0/e2*hyin(m,1:index-1);
        exin(m,index:end)=beta(m)/k0*eta0/eBar*hyin(m,index:end);
    end
    
    for nIdx=1:2*N+1
        hyout(nIdx,:)=exp(-1i*kxn(nIdx).*x);
        exout(nIdx,:)=gamma(nIdx)/k0*eta0/e1*hyout(nIdx,:);
    end
            
    hin=hyin;
    ein=exin;
    hout=hyout;
    eout=exout;
    
%%%% TE polarization
else
    hxin=zeros(M,length(x));
    eyin=zeros(M,length(x));
    hxout=zeros(2*N+1,length(x));
    eyout=zeros(2*N+1,length(x));

    for m=1:1:M       
        hxin(m,1:index-1)=A(m)*cos(ka(m)*(x(1:index-1)-a/2))+B(m)*sin(ka(m)*(x(1:index-1)-a/2));
        hxin(m,index:end)=C(m)*cos(ks(m)*(x(index:end)-(a+period)/2))+D(m)*sin(ks(m)*(x(index:end)-(a+period)/2));        
        eyin(m,1:index-1)=-(k0/beta(m))*eta0*u2*hxin(m,1:index-1);
        eyin(m,index:end)=-(k0/beta(m))*eta0*uBar*hxin(m,index:end);
    end
    
    for nIdx=1:2*N+1
        hxout(nIdx,:)=exp(-1i*kxn(nIdx).*x);
        eyout(nIdx,:)=-(k0/gamma(nIdx))*eta0*u1*hxout(nIdx,:);
    end
    
    hin=hxin;
    ein=eyin;
    hout=hxout;
    eout=eyout;
end

%% Calculate matrixes, reflection and transmission coefficient
[powRefl, powTrans, ~, ~, ~] = solveRT_OA(H, E, beta, gamma, tg, lambda, e1, u1, a+s, theta, N, M, isTM);

% matrix: phi
phi = diag(exp(-1i*beta*tg), 0);	% phi is a diagonal matrix
% matrix: rho
rho=(((eye(M)+(H\eye(M))*E))\eye(M))*(eye(M)-(H\eye(M))*E);
% matrix: Zin, R
Zin=E*(eye(M)+phi*rho*phi)*((eye(M)-phi*rho*phi)\eye(M))*(H\eye(M));
R=((Zin+eye(2*N+1))\eye(M))*(Zin-eye(2*N+1));
% vector: r
r=(R*deltaN0')';
% vector: aa, tou
aa=2*phi*(((Zin\eye(M)+eye(2*N+1))*E*(eye(M)+phi*rho*phi))\eye(M))*deltaN0';
tau=E*(eye(M)+rho)*aa;
% matrix: T
T=2*E*(eye(M)+rho)*phi*(((Zin\eye(M)+eye(2*N+1))*E*(eye(M)+phi*rho*phi))\eye(M));  % Good for 2*N+1=M

%% Calculate field profile
aarho=rho*aa;
% construction of E, H field in region II
index1=find(z<-tg);
index2=find(z<0);
z1=z(1:index1(end));
z2=z(index1(end)+1:index2(end));
z3=z(index2(end)+1:end);

Htemp=zeros(length(z2),M);
Etemp=zeros(length(z2),M);
for m=1:M
    Htemp(:,m)=aa(m)*exp(-1i*beta(m)*z2)-aarho(m)*exp(1i*beta(m)*z2);
    Etemp(:,m)=aa(m)*exp(-1i*beta(m)*z2)+aarho(m)*exp(1i*beta(m)*z2);
end
H2=Htemp*hin;
E2=Etemp*ein;

% construction of E, H field in region I
Htemp=zeros(length(z1),2*N+1);
Etemp=zeros(length(z1),2*N+1);
for n=1:2*N+1
    Htemp(:,n)=exp(1i*gamma(n)*(z1+tg))*conj(r(n));  
    Etemp(:,n)=exp(1i*gamma(n)*(z1+tg))*conj(r(n)); 
end
HincidentTempz=exp(1i*kz0*(z1+tg));
HincidentTempx=exp(-1i*kx0*x);                       
Hincident=HincidentTempz'*HincidentTempx;

H1=Hincident-Htemp*hout;
if isTM                        
    E1=eta0/e1*(gamma((1+end)/2))/k0*Hincident+Etemp*eout;
else
    E1=-eta0*u1*k0/(gamma((1+end)/2))*Hincident+Etemp*eout;
end

% construction of E, H field in region III
Htemp=zeros(length(z3),2*N+1);
Etemp=zeros(length(z3),2*N+1);
for n=1:2*N+1
    Htemp(:,n)=exp(-1i*gamma(n)*z3)*tau(n);
    Etemp(:,n)=exp(-1i*gamma(n)*z3)*tau(n);
end
H3=Htemp*hout;
E3=Etemp*eout;

H_all=[H1; H2; H3];
E_all=[E1; E2; E3];

if isTM
    incidentWaveAmp=[1 eta0/e1*gamma((1+end)/2)/k0];
else
    incidentWaveAmp=[1 -eta0*u1*k0/(gamma((1+end)/2))];
end

end

