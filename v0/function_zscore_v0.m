%Medical Physics Department at Bariloche Atomic Center, Argentina.
%Author: Damian Dellavale (dellavale@cab.cnea.gov.ar, dellavaledamian@gmail.com).
%Project: Window functions.
%Date: 14/12/2015.

%Description:
%In this script the input signal is z-score normalized in order to have
%zero mean and unit variance.

%Tree of dependencies:
%None.

%Reference:
%Velarde O, Urdapilleta E, Mato G, and Dellavale D (2019), Bifurcation
%structure determines different phase-amplitude coupling patterns in the
%activity of biologically plausible neural networks, NeuroImage, In Press,
%(DOI: ...)

function [normSignal] = function_zscore_v0(rawSignal)
%==========================================================================
%Inputs:
%rawSignal -> Raw signal to be normalized (column array: samples x channels).

%Outputs:
%normSignal -> Normalized signal (column array: samples x channels).
%==========================================================================

%Argument completion ------------------------------------------------------
if (nargin < 1), error('MATLAB:function_zscore','Input argument error.'); end
%--------------------------------------------------------------------------

%Check the input arguments ------------------------------------------------
%assert(max(size(rawSignal))==size(rawSignal,1), 'Input argument error in function "function_zscore": The signal must be a column array.');

if (max(size(rawSignal))>size(rawSignal,1)), 
    warning('MATLAB:function_zscore',...
            ['The normalization is performed across the rows. '...
             'The input matrix has Nrows < Ncol. '...
             'Consider traspose the input matrix to implement the normalization across the columns.']);
end
%--------------------------------------------------------------------------

%Default values of the outputs --------------------------------------------
%--------------------------------------------------------------------------

%Parameters ---------------------------------------------------------------
%Set the dimention along to which compute the z-score normalization.
dim = 1;
%--------------------------------------------------------------------------

%Normalization in order to have zero mean and unit variance ---------------

%Zero mean normalization.
offset = mean(rawSignal,dim);
%offset = repmat(offset,size(rawSignal,dim),1);
offset = offset(ones(size(rawSignal,dim),1),:); %This one is faster.
normSignal = rawSignal - offset;

%Unit variance normalization.
%Compute the standard deviation normalized with "N-1", where "N" is the number of samples. 
flag = 0;
stdDev  = std(normSignal,flag,dim);
%stdDev = repmat(stdDev,size(normSignal,dim),1);
stdDev = stdDev(ones(size(normSignal,dim),1),:); %This one is faster.
normSignal = normSignal./stdDev;

%Verify the normalizations:
%mean(normSignal,dim)
%var(normSignal,dim)
%--------------------------------------------------------------------------

end %function

