%Medical Physics Department at Bariloche Atomic Center, Argentina.
%Author: Damian Dellavale (dellavale@cab.cnea.gov.ar, dellavaledamian@gmail.com).
%Project: Time Locked Index (TLI).
%Date: 10/06/2019.

%Description:
%In this script the Time Locked Index (TLI) is computed.
%The zero-crossings of the LF-phase time series are used
%to find the indices for the peaks of LF and HF time series.

%Tree of dependencies:
%function_TimeLockedIndex_v0.m
% function_FDF_v0.m
%  function_window_v0.m
% function_zscore_v0.m

%Changes from previous versions:
%None.

%Reference:
%Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation
%structure determines different phase-amplitude coupling patterns in the
%activity of biologically plausible neural networks, NeuroImage, In Press,
%(DOI: ...)

function [TLIout] = function_TimeLockedIndex_v0(LFsignal, HFsignal, TLIcfg)
%==========================================================================
%Inputs:
%LFsignal  -> Band-pass filtered signal using the filter configuration TLIcfg.BPFcfg_LF [signal] (column array: samples x 1). 
%HFsignal  -> Band-pass filtered signal using the filter configuration TLIcfg.BPFcfg_HF [signal] (column array: samples x 1). 
%TLIcfg -> Time Locked Index configuration (structure array).
%          'BPFcfg' -> Band-Pass Filter configuration (structure array).
%          'abs_HFSignaltimeLockedHFpeaks' -> Flag for using the absolute value in detecting the peaks of the HF signal and,
%                                             to compute the HF signal timed-locked averaged to the FAST oscillations peaks 
%                                             "HFSignaltimeLockedHFpeaks".
%          'abs_HFSignaltimeLockedLFpeaks' -> Flag for using the absolute value to compute the HF signal timed-locked averaged
%                                             to the SLOW oscillations peaks "HFSignaltimeLockedLFpeaks".
%          'LFphase' -> Phase of the LF ("peaks" or "troughs") to define the time interval interval to find the HF peaks.
%          'NT'      -> Number of LF periods around the fast oscillation peaks [none] (scalar).
%          'plot'    -> Flag to plot the signals.
%                       0 => Do not plot.
%                       otherwise => Plot the signals involved in the TLI computation.
%          'fs'      -> Sampling rate [Hz] (scalar).
%
%Outputs:
%TLIout -> Outut structure (structure).
%HFSignaltimeLockedHFpeaks  -> HF signal timed-locked averaged to the FAST oscillations peaks [signal] (column array: samples x 1).
%HFSignaltimeLockedLFpeaks  -> HF signal timed-locked averaged to the SLOW oscillations peaks [signal] (column array: samples x 1).
%TLI                       -> Time-Locked Index (1 x 1).
%indPeak_HF                -> Peak indices corresponding to the HF fast oscillations.
%indPeak_LF                -> Peak indices corresponding to the LF slow oscillations.
%sampleT                   -> Epoch length (multiple of the LF signal period).
%==========================================================================

%Argument completion ------------------------------------------------------
if (nargin < 3)||isempty(LFsignal)||isempty(HFsignal)||isempty(TLIcfg),...
   error('MATLAB:function_TimeLockedIndex','Input argument error.');
end

if isfield(TLIcfg.BPFcfg_LF, 'f1') && isfield(TLIcfg.BPFcfg_LF, 'f2'),
    
    %Compute the center frequency.
    TLIcfg.BPFcfg_LF.f0 = ( TLIcfg.BPFcfg_LF.f2 + TLIcfg.BPFcfg_LF.f1 ) / 2; %Arithmetic mean. 
    %TLIcfg.BPFcfg_LF.f0 = sqrt( TLIcfg.BPFcfg_LF.f2 * TLIcfg.BPFcfg_LF.f1 ); %Geometric mean.
    %Ref: https://en.wikipedia.org/wiki/Center_frequency
    
    %Compute the bandwidth.
    TLIcfg.BPFcfg_LF.Bw = TLIcfg.BPFcfg_LF.f2 - TLIcfg.BPFcfg_LF.f1;
    
elseif ~isfield(TLIcfg.BPFcfg_LF, 'f0') || ~isfield(TLIcfg.BPFcfg_LF, 'Bw'),
    
    error('MATLAB:function_TimeLockedIndex','Error in the BPF configuration (BPFcfg_LF).');
    
end

if isfield(TLIcfg.BPFcfg_HF, 'f1') && isfield(TLIcfg.BPFcfg_HF, 'f2'),
    
    %Compute the center frequency.
    TLIcfg.BPFcfg_HF.f0 = ( TLIcfg.BPFcfg_HF.f2 + TLIcfg.BPFcfg_HF.f1 ) / 2; %Arithmetic mean. 
    %TLIcfg.BPFcfg_HF.f0 = sqrt( TLIcfg.BPFcfg_HF.f2 * TLIcfg.BPFcfg_HF.f1 ); %Geometric mean.
    %Ref: https://en.wikipedia.org/wiki/Center_frequency
    
    %Compute the bandwidth.
    TLIcfg.BPFcfg_HF.Bw = TLIcfg.BPFcfg_HF.f2 - TLIcfg.BPFcfg_HF.f1;
    
elseif ~isfield(TLIcfg.BPFcfg_HF, 'f0') || ~isfield(TLIcfg.BPFcfg_HF, 'Bw'),
    
    error('MATLAB:function_TimeLockedIndex','Error in the BPF configuration (BPFcfg_HF).');
    
end
%--------------------------------------------------------------------------

%Check the input arguments ------------------------------------------------
assert(size(LFsignal,2)==1, 'Input argument error in function "function_TimeLockedIndex": The LFsignal must be a column array.');
assert(size(HFsignal,2)==1, 'Input argument error in function "function_TimeLockedIndex": The HFsignal must be a column array.');

assert(isstruct(TLIcfg), 'Input argument error in function "function_TimeLockedIndex": TLIcfg must be a structure array.');

if ~isfield(TLIcfg, 'NT')||isempty(TLIcfg.NT), TLIcfg.NT = 1; end %Default value.
if ~isfield(TLIcfg, 'fs')||isempty(TLIcfg.fs), TLIcfg.fs = 1; end %Default value.
%--------------------------------------------------------------------------

%Default values of the outputs --------------------------------------------
HFSignaltimeLockedHFpeaks = NaN;
HFSignaltimeLockedLFpeaks = NaN;
TLI = NaN;
%--------------------------------------------------------------------------

%Parameters ---------------------------------------------------------------
%Compute the number of samples for the lower cutoff frequency.
lowestFrec = TLIcfg.BPFcfg_LF.f0 - TLIcfg.BPFcfg_LF.Bw/2;
sampleT = round( TLIcfg.NT * TLIcfg.fs / lowestFrec );
halfSampleT = round(sampleT/2);

%Compute the length of the raw signal.
Nsamples = size(LFsignal,1);
%--------------------------------------------------------------------------

%Z-score normalization ----------------------------------------------------  
%in order to have zero mean and unit variance (unit amplitude) in HF and LF signals.

HFsignal = function_zscore_v0(HFsignal);
LFsignal = function_zscore_v0(LFsignal);

%Verify the normalizations:
%mean(HFsignal,1)
%mean(LFsignal,1)
%var(HFsignal,1)
%var(LFsignal,1)
%--------------------------------------------------------------------------

%Find the index for the peak corresponding to phase = 0 or pi rad. --------

%Compute the LF-phase time series.
LFphase = angle(hilbert(LFsignal)); %[rad] range: [-pi,pi]
%NOTE: The LFphase time series is only used to identify the peaks or troughs
%of the LF signal. As a consequence, here we can disregard the transient 
%response due to the Hilbert transform.
%The Hilbert transform is computed in the frequency domain in the Matlab
%script hilbert.m via "fft" and "ifft" built-in functions. 

%Find zero-crossings in LF-phase time series.
s1 = LFphase(1:Nsamples-1);
s2 = LFphase(2:Nsamples);

switch lower(TLIcfg.LFphase), %Switch for LFphase.
    case 'peaks',
        pulseZC = (s1.*s2)<0 & abs(s1-s2)<pi; %Peaks: LF phase = 0 rad,
        %i.e., the maximum (absolute) amplitude of the LF signal in each period.
    case 'troughs',
        pulseZC = (s1.*s2)<0 & abs(s1-s2)>=pi; %Troughs: LF phase = pi rad,
        %i.e., the minimum (absolute) amplitude of the LF signal in each period.
    otherwise,
        error('MATLAB:function_TimeLockedIndex',...
              ['The option: ',TLIcfg.LFphase,' does not exist.']);
end %Switch for LFphase.

%Re-arrange the pulses in order to emulate a causal zero-crossing detection.
dim = 1;
pulseZC = cat(dim, 0, pulseZC);

%Compute the indices for the peaks of the LF signal (zeroes of the LF-phase signal).
indPeak_LF = find(pulseZC > 0);

% %DEBUGGING: Show the signals ---
% function_showData_v1([LFphase/max(abs(LFphase)), LFsignal/max(abs(LFsignal))]);
% 
% figure, plot(LFphase,'.-b')
% hold on, plot(pulseZC,'or')
% 
% figure, plot(LFsignal,'.-b')
% hold on, plot(indPeak_LF,LFsignal(indPeak_LF),'+r','Markersize',20)
% %---

%Compute the number of LF peaks.
Npeak = length(indPeak_LF);
%--------------------------------------------------------------------------

%Compute the auxiliary signals according to the flags for the absolute value.
if TLIcfg.abs_HFSignaltimeLockedHFpeaks,
    HFsignal_aux1 = abs(HFsignal);        
else
    HFsignal_aux1 = HFsignal;        
end

if TLIcfg.abs_HFSignaltimeLockedLFpeaks,
    HFsignal_aux2 = abs(HFsignal);        
else
    HFsignal_aux2 = HFsignal;        
end
%--------------------------------------------------------------------------

%Compute the peaks of HF --------------------------------------------------

%Memory pre-allocation for speed up the loop.
indPeak_HF = NaN(Npeak-1,1);
for ii=1:+1:Npeak-1, %Loop across all the "Npeaks-1" intervals between the LF peaks.
    
    %Find the index for the peak corresponding to the maximum (absolute)
    %amplitude of the HF signal in each period of the LF signal.
    [~, indPeak_HF(ii)] = max( HFsignal_aux1(indPeak_LF(ii):indPeak_LF(ii+1)) );

    indPeak_HF(ii) = indPeak_HF(ii) + indPeak_LF(ii) - 1;
    
end %Loop across the peaks.

%Note that it results:
%length(indPeak_HF) = length(indPeak_LF) - 1 = Npeak - 1;
%--------------------------------------------------------------------------
    
%Extract the epochs centered at the time points corresponding to the
%FAST and SLOW OSCILLATION PEAKS ------------------------------------------

%Compute the valid indices for the HF and LF peaks.
JJstart = NaN;
JJend = NaN;

for jj=1:+1:Npeak-1, %Loop across all the HF peaks: length(indPeak_HF) = length(indPeak_LF) - 1 = Npeak - 1;
    if (indPeak_LF(jj)-halfSampleT > 0) && (indPeak_HF(jj)-halfSampleT > 0),
        JJstart = jj;
        break,
    end
end
for jj=Npeak-1:-1:1, %Loop across all the HF peaks: length(indPeak_HF) = length(indPeak_LF) - 1 = Npeak - 1;
    if (indPeak_LF(jj)+halfSampleT < Nsamples) && (indPeak_HF(jj)+halfSampleT < Nsamples),
        JJend = jj;
        break,
    end
end
if isnan(JJstart)||isnan(JJend),
    warning('MATLAB:function_TimeLockedIndex',...
            'Undefined bounds for the indices of the HF/LF peaks.');
    return,    
end

%Memory pre-allocation for speed up the loop.
HFSignaltimeLockedHFpeaks = zeros(2*halfSampleT+1,1);
HFSignaltimeLockedLFpeaks = HFSignaltimeLockedHFpeaks;
for jj=JJstart:+1:JJend, %Loop over the peak indices.
    
    HFSignaltimeLockedHFpeaks = HFSignaltimeLockedHFpeaks + HFsignal_aux1(indPeak_HF(jj)-halfSampleT:indPeak_HF(jj)+halfSampleT);
    HFSignaltimeLockedLFpeaks = HFSignaltimeLockedLFpeaks + HFsignal_aux2(indPeak_LF(jj)-halfSampleT:indPeak_LF(jj)+halfSampleT);

end
    
%Compute the average of the time-locked epochs.
HFSignaltimeLockedHFpeaks = HFSignaltimeLockedHFpeaks / (JJend-JJstart+1);
HFSignaltimeLockedLFpeaks = HFSignaltimeLockedLFpeaks / (JJend-JJstart+1);
%--------------------------------------------------------------------------

%Compute the Time-Locked Index --------------------------------------------
TLI = ( max(HFSignaltimeLockedLFpeaks) - min(HFSignaltimeLockedLFpeaks) ) /...
      ( max(HFSignaltimeLockedHFpeaks) - min(HFSignaltimeLockedHFpeaks) );
%
%If absolute values are used, apply the corresponding scale factor so the
%resulting TLI is constrained to the interval [0,1].
if TLIcfg.abs_HFSignaltimeLockedHFpeaks, TLI = TLI/2; end
if TLIcfg.abs_HFSignaltimeLockedLFpeaks, TLI = 2*TLI; end
%--------------------------------------------------------------------------

%Synthesize the output structure ------------------------------------------
TLIout = struct(...
'HFSignaltimeLockedHFpeaks', HFSignaltimeLockedHFpeaks,...
'HFSignaltimeLockedLFpeaks', HFSignaltimeLockedLFpeaks,...
'TLI', TLI,... 
'indPeak_HF', indPeak_HF,... 
'indPeak_LF', indPeak_LF,... 
'sampleT', sampleT);
%--------------------------------------------------------------------------

end %function

