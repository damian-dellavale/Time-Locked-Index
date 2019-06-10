%Medical Physics Department at Bariloche Atomic Center, Argentina.
%Author: Damian Dellavale.
%Project: Window functions.
%Date: 01/01/2016.

%Description:
%In this script we synthesize the weighting function for windowing in time and frequency domain.

%Tree of dependencies:
%none.

%Reference:
%Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation
%structure determines different phase-amplitude coupling patterns in the
%activity of biologically plausible neural networks, NeuroImage, In Press,
%(DOI: ...)

function windowFunction = function_window_v0(windowParam, windowLength)
%==========================================================================
%Inputs:
%windowParam    -> Parameters of the window function (structure array).
%                  name  -> Name defining the window type (string).
%                  alpha -> Parameter for the gausswin: alpha is inversely proportional to the standard
%                           deviation of a Gaussian random variable.
%                  sflag -> Sampling parameter for hamming and hann windows (string).
%                           "sflag" can be either "periodic" or "symmetric" (the default).
%                           The "periodic" flag is useful for DFT/FFT purposes, such as in spectral analysis.
%                           The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal windowed
%                           with a periodic window to have perfect periodic extension. When "periodic" is specified,
%                           hamming/hann computes a length "windowLength+1" window and returns the first "windowLength" points.
%                           When using windows for filter design, the "symmetric" flag should be used. 
%                  r     -> A Tukey window is a rectangular window with the first and last 100*r/2 percent of the samples equal
%                           to parts of a cosine. r is a real number between 0 and 1.
%                           If you input r <= 0, you obtain a rectwin window.
%                           If you input r >= 1, you obtain a hann window.
%windowLength   -> Length of the window function.
%
%Outputs:
%windowFunction -> Synthesized window function.
%==========================================================================

%Argument completion ------------------------------------------------------
if (nargin < 2)||isempty(windowParam)||isempty(windowLength),...
   error('MATLAB:function_window','Input argument error.');
end
%--------------------------------------------------------------------------

%Check the input arguments ------------------------------------------------
assert(isstruct(windowParam), 'Input argument error in function "function_window": windowParam must be a structure array.');
%--------------------------------------------------------------------------

switch lower(windowParam.name), %Switch for window's selection.
    case 'gausswin', %Gaussian window:
        
        windowFunction = gausswin(windowLength, windowParam.alpha);
        %windowFunction = window(@gausswin, windowLength, windowParam.alpha);
        
%         warning('MATLAB:function_window', ['The Gaussian window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']);

    case 'hamming', %Hamming window:
        
        windowFunction = hamming(windowLength, windowParam.sflag);
        %windowFunction = window(@hamming, windowLength, windowParam.sflag);
        
%         warning('MATLAB:function_window', ['The Hamming window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']);
                                       
    case 'hann', %Hann (Hanning) window:
        
        windowFunction = hann(windowLength, windowParam.sflag);
        %windowFunction = window(@hann, windowLength, windowParam.sflag);
        
%         warning('MATLAB:function_window', ['The Hann window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']);   

    case 'tukey', %Tukey (tapered cosine) window:
        
        windowFunction = tukeywin(windowLength, windowParam.r);
        %windowFunction = window(@tukeywin, windowLength, windowParam.r);
        
%         warning('MATLAB:function_window', ['The Tukey window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']); 

    otherwise, %Rectangular window:
        %This function is provided for completeness; a rectangular window is equivalent to no window at all.
        windowFunction = rectwin(windowLength);
        %windowFunction = ones(windowLength,1);
end %Switch for window's selection.

end %function
