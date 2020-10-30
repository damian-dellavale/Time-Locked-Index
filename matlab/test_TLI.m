%Medical Physics Department at Bariloche Atomic Center, Argentina.
%Authors: Osvaldo Velarde (osva.m.velarde@gmail.com), 
%         Damian Dellavale (dellavale@cab.cnea.gov.ar, dellavaledamian@gmail.com).
%Project: Test script for the Time Locked Index (TLI) using the Van der Pol oscillator.
%Date: 10/06/2019.

%Description:
%The TLI metric is computed for the dynamics of the Van der Pol oscillator.
%
%Van der Pol differential equation:
%d^2 (x) - mu (1-x^2) d(x) + w^2 * x = 0

%Reference:
%Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation
%structure determines different phase-amplitude coupling patterns in the
%activity of biologically plausible neural networks, NeuroImage, In Press,
%(DOI: ...)

%%
clear global        %Clear all global variables.
clear all           %Clear all local variables.
close all force     %Close all figures (including the wvtool figures).
format long         %Show numbers in long format.
clc                 %Clear screen.
restoredefaultpath
%pathtool

%% Access path to the functions.

pathFunctions1 = './';

addpath(pathFunctions1,...
        '-begin');

%% Parameters of the Van der Pol equation.

%Fundamental frequency of the oscillator at mu=0.
f = 10; %[Hz]

%Angular frequency.
w = 2*pi*f; %[rad/sec.] 

%Non linear parameter of the oscillator.
mu = (0:+2:250);

%% Parameters for the solver.

%Time resolution.
dt = 0.001; %[sec.]

%Sampling frequency.
fs = 1/dt; %[Hz]

%Time lenght of the simulation.
TIME_SIM = 6; %[sec.]

%Define the transient of the solver: one period of the fundamental frequency of the oscillator.
SOLVER_TRANSIENT = fs/f;    
    
%% Frequency bands of interest.

FREQ_BANDS = [1, 15;...   %Low-frequency band [Hz]: It includes the fundamental 
                          %frequency of the Van der Pol oscillator for all 
                          %the explored values of the parameter mu.
                          %
              20, 100;... %High-frequency band [Hz]: It comprises several 
                          %harmonics of the fundamental frequency.
              ];    

%% Filtering parameters.
%NOTE: Change this section to use your preferred filtering function.

%Define the band-pass filter to compute the low-frequency time series.
indLF = 1;
f_LF = (FREQ_BANDS(indLF,2)+FREQ_BANDS(indLF,1))/2; %[Hz] Center frequency.
Bw_LF = FREQ_BANDS(indLF,2)-FREQ_BANDS(indLF,1); %[Hz] Null-to-null bandwidth.
%Configuration structure for filtering in the low-frequency band.
BPFcfg_LF = struct('f0',f_LF,...
                   'Bw',Bw_LF,...
                   'zeropadding',0,...
                   'freqWindowParam',struct('name','hann','sflag','symmetric'),...
                   'timeWindowParam',struct('name','hann','sflag','symmetric'),...
                   'conv','circular',...
                   'causal',false,...
                   'fs',fs);

%Define the band-pass filter to compute the high-frequency time series.
indHF = 2;
f_HF = (FREQ_BANDS(indHF,2)+FREQ_BANDS(indHF,1))/2; %[Hz] Center frequency.
Bw_HF = FREQ_BANDS(indHF,2)-FREQ_BANDS(indHF,1); %[Hz] Null-to-null bandwidth.
%Configuration structure for filtering in the high-frequency band.
BPFcfg_HF = struct('f0',f_HF,...
                   'Bw',Bw_HF,...
                   'zeropadding',0,...
                   'freqWindowParam',struct('name','hann','sflag','symmetric'),...
                   'timeWindowParam',struct('name','hann','sflag','symmetric'),...
                   'conv','circular',...
                   'causal',false,...
                   'fs',fs);

%% Input parameters to compute the Time Locked Index.

TLIcfg = struct('BPFcfg_LF',BPFcfg_LF,...
                'BPFcfg_HF',BPFcfg_HF,...
                'abs_HFSignaltimeLockedHFpeaks',false,...
                'abs_HFSignaltimeLockedLFpeaks',false,...
                'LFphase','peaks',...
                'NT',1,...
                'plot',false,...
                'fs',fs);
            
%Define the noise level.
NOISE_LEVEL = 0.1; %Fraction of the signal's standard deviation.

%Control random number generation.
%seed = 0;
%rng(seed);
%Ref: https://www.mathworks.com/help/matlab/ref/rng.html

%% Implement the simulations for the Van der Pol oscillator.

signal = function_vanDerPolSolver_v0(mu, w, TIME_SIM, dt);
signal(1:SOLVER_TRANSIENT,:) = []; %Remove the transient of the solver.

Ns = size(signal,1);   %Compute the number of samples of the input signal.
t = (0:+1:Ns-1).' / fs; %[sec] %Synthesize the time vector.

%% Compute the TLI as a function of the parameter mu.

%Memory pre-allocation to speed up the loop.
TLI = NaN(1,length(mu));
for ii=1:+1:length(mu), %Loop across the mu values.
    
    display(['Processing: mu = ',num2str(mu(ii)),'; ',num2str(ii),' out of ',num2str(length(mu)),'.'])
    
    %Load the raw signal for the current mu value.
    rawSignal = signal(:,ii);

    %----------------------------------------------------------------------
    %Additive white Gaussian noise.

    %Compute the standard deviation, normalized with "N-1", where "N" is the number of samples. 
    flag = 0;
    stdDevSignal  = std(rawSignal,flag,1);

    %Include additive white Gaussian noise.
    noise = (NOISE_LEVEL * stdDevSignal) * rand(size(rawSignal));
    rawSignal = rawSignal + noise;
    %----------------------------------------------------------------------

    %Z-score normalization of the input signal.
    rawSignal = function_zscore_v0(rawSignal);
           
    %Reflect the time series to minimize edge artifacts due to the transient response of the BPFs.
    %Ref: Mike X. Cohen, Analyzing Neural Time Series Data, Theory and Practice, MIT press, 2014.
    %      Figure 7.3, p 78
    %      14.9 Filtering Each Trial versus Filtering Concatenated Trials, p 190 
    %      Figure 14.11, p 191
    %
    rawSignal = [rawSignal(end:-1:1); rawSignal; rawSignal(end:-1:1)];

    %Compute the Band-Pass Filtering.
    %NOTE: Change this piece of code to use your preferred filtering function.
    [FDFout] = function_FDF_v0(rawSignal, BPFcfg_LF);
    LFSIGNAL = FDFout.filteredSignal;
    [FDFout] = function_FDF_v0(rawSignal, BPFcfg_HF);
    HFSIGNAL = FDFout.filteredSignal;        
    
    %Restore the length of the signals to remove the transient of the filtering.
    rawSignal = rawSignal(Ns:end-Ns-1,:);
    LFSIGNAL = LFSIGNAL(Ns:end-Ns-1,:);
    HFSIGNAL = HFSIGNAL(Ns:end-Ns-1,:);

    %Compute the TLI.
    [TLIout] = function_TimeLockedIndex_v0(LFSIGNAL, HFSIGNAL, TLIcfg);
    TLI(ii) = TLIout.TLI;
    
    if ii==1 || ii==length(mu),
        figure, hold on, 
        plot(t,rawSignal,'-k')
        plot(t,LFSIGNAL,'-b')
        plot(t,HFSIGNAL,'-r')
        legend({'Raw signal','LF signal','HF signal'},'location','best','interpreter','none')
        box on, grid off,
        axis tight,
        set(gca,'fontsize',24)
        xlabel(['Time [sec.]'],'interpreter','none')
        ylabel(['Amplitude [arb. units]'],'interpreter','none')
        title(['$\mu / \omega$ = ',num2str(mu(ii)/w,'%1.3g')],'interpreter','LaTex')
    end

end %Loop across the mu values.

%Plot the TLI values as a function of the parameter mu.
figure, plot(mu/w,TLI,':ob','Markersize',10)
box on, grid on,
axis tight, ylim([0, 1.1]);
set(gca,'fontsize',24)
xlabel('$\mu / \omega$ [arb. units]','interpreter','LaTex')
ylabel('$TLI$ [arb. units]','interpreter','LaTex')
