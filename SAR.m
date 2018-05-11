%%% Author: Aditya Jain %%%%
%%% Date: 6th April, 2018 %%%
%%% Topic: Synthetic Aperture Radar (SAR) %%%%

clear
close all
clc
tic
%% Variable Declaration

bw = 20e+6;  % Bandwidth in hz 
c = 3e+8;    % Speed of light
Tf = 34e-6;  % Fast time in seconds
k = bw/Tf;   % chirp factor
vp = 6600;   % Speed of platform in m/sec
fc = 1.257e+6; % carrier frequency
R = 854000;  % Slant range
Nslow = 1024; % No of slow stops
lambda = c/fc; % Wavelength

%% 
Fmax = (2*(vp^2)*Tf)/(lambda*R);
Fs = 3*bw;   % sampling frequency
delt = 1/Fs; % resolution on time axis

% time axis'
tfast = 0:delt:Tf-delt;  % fast time axis
rangeaxis = (-Tf+delt:delt:Tf-delt)*(c/2);     % range axis
Nfast = Tf/delt;   % no of points on the fast axis

tslow = 0:Tf:(Nslow-1)*Tf; % slow time axis
craxis = (-(Nslow-1)*Tf:Tf:(Nslow-1)*Tf)*vp;   % cross-range axis

Atgt = 1;         % magnitude of the received signal
RTarget = 2000;   % target distance in m 
crTarget = 100;   % cross-range target distance

SrxDelayed = zeros(Nslow, Nfast);     % Shifted received signal
SrxZeroShift = zeros(Nslow, Nfast);   % Received signal without any shift

for m = 1:Nslow
    for n = 1:Nfast
        
        t2 = tslow(m);
        t1 = tfast(n);
        
        % After down-conversion
        SrxDelayed(m,n) = Atgt*exp(1i*pi*k*(t1 - 2*RTarget/c)^2)*exp((-1i*2*pi*(vp)^2*(t2 - crTarget/vp)^2)/(lambda*R));
        SrxZeroShift(m,n) = Atgt*exp(1i*pi*k*t1^2)*exp((-1i*2*pi*vp^2*t2^2)/(lambda*R));
        
    end    
end


% Applying 2D compression filter by using 2D cross-correlation
SrxPC = xcorr2(SrxDelayed, SrxZeroShift);

% Plotting
%%
SrxPC_db = 20*log10(abs(SrxPC));
cmax = max(max(SrxPC_db));
imagesc(rangeaxis, craxis, SrxPC_db,[cmax-30 cmax]);
% imagesc(rangeaxis, craxis, SrxPC_db);
xlabel('Range (in m)')
ylabel('Cross-Range (in m)')
title('Synthetic Aperture Radar (SAR)')

toc