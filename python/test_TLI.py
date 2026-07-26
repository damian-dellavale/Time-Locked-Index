import numpy as np
import scipy.signal.windows as windows
from scipy.signal import hilbert
from scipy.integrate import solve_ivp
from scipy.fft import fft, ifft, fftshift
import matplotlib.pyplot as plt
import math

# ==========================================
# 1. Z-Score Normalization
# ==========================================
def function_zscore_v0(rawSignal):
    """
    Normalizes the input signal to have zero mean and unit variance.
    """
    mean_val = np.mean(rawSignal, axis=0)
    std_val = np.std(rawSignal, axis=0, ddof=1)
    # Prevent division by zero if std is zero
    std_val[std_val == 0] = 1.0
    return (rawSignal - mean_val) / std_val


# ==========================================
# 2. Window Functions
# ==========================================
def function_window_v0(windowParam, windowLength):
    """
    Synthesizes the weighting function for windowing in time and frequency domains.
    """
    windowLength = int(windowLength)
    if windowLength <= 0:
        return np.array([])

    name = windowParam.get('name', 'rectwin').lower()
    
    if name == 'gausswin':
        alpha = windowParam.get('alpha', 2.5)
        std = (windowLength - 1) / (2 * alpha) if windowLength > 1 else 1e-10
        return windows.gaussian(windowLength, std=std)
    
    elif name == 'hamming':
        sflag = windowParam.get('sflag', 'symmetric')
        sym = (sflag.lower() == 'symmetric')
        return windows.hamming(windowLength, sym=sym)
        
    elif name == 'hann':
        sflag = windowParam.get('sflag', 'symmetric')
        sym = (sflag.lower() == 'symmetric')
        return windows.hann(windowLength, sym=sym)
        
    elif name == 'tukey':
        r = windowParam.get('r', 0.5)
        return windows.tukey(windowLength, alpha=r)
        
    else:
        return np.ones(windowLength)


# ==========================================
# 3. Frequency Domain Filtering (FDF)
# ==========================================
def function_FDF_v0(signal, FDFcfg):
    """
    Implements Frequency Domain Filtering (FDF).
    """
    if signal.ndim == 1:
        signal = signal.reshape(-1, 1)
        
    fs = FDFcfg.get('fs', 1.0)
    zeropadding = FDFcfg.get('zeropadding', 0)
    conv_type = FDFcfg.get('conv', 'circular').lower()
    
    # Determine cutoff frequencies
    if 'f0' in FDFcfg and 'Bw' in FDFcfg:
        f1 = FDFcfg['f0'] - FDFcfg['Bw'] / 2.0
        f2 = FDFcfg['f0'] + FDFcfg['Bw'] / 2.0
    else:
        f1 = FDFcfg.get('f1', 0.0)
        f2 = FDFcfg.get('f2', fs / 2.0)
    
    Ns, Nch = signal.shape
    
    # Calculate next power of 2 length
    next_pow2 = math.ceil(math.log2(Ns)) if Ns > 0 else 0
    nfft = int(2 ** (next_pow2 + zeropadding))
    
    if zeropadding < 0:
        signal_padded = signal[:nfft, :]
    else:
        signal_padded = np.vstack([signal, np.zeros((nfft - Ns, Nch))])
        
    freqs = np.fft.fftfreq(nfft, d=1/fs)
    
    # Construct frequency mask
    H = np.zeros(nfft)
    pass_idx = (np.abs(freqs) >= f1) & (np.abs(freqs) <= f2)
    H[pass_idx] = 1.0 
    
    FFT_sig = fft(signal_padded, n=nfft, axis=0)
    HH = H[:, np.newaxis]
    filtered_FFT = FFT_sig * HH
    
    filteredSignal = np.real(ifft(filtered_FFT, n=nfft, axis=0))
    
    # Restore signal length if circular convolution
    if conv_type == 'circular' and (nfft >= Ns):
        filteredSignal = filteredSignal[:Ns, :]
        
    return {'filteredSignal': filteredSignal}


# ==========================================
# 4. Time Locked Index (TLI)
# ==========================================
def function_TimeLockedIndex_v0(LFsignal, HFsignal, TLIcfg):
    """
    Computes the Time Locked Index (TLI) using zero-crossings of the LF-phase.
    """
    fs = TLIcfg.get('fs', 1.0)
    NT = TLIcfg.get('NT', 1.0)
    LFphase_type = TLIcfg.get('LFphase', 'peaks').lower()
    
    f0_LF = (TLIcfg['BPFcfg_LF']['f1'] + TLIcfg['BPFcfg_LF']['f2']) / 2.0
    Bw_LF = TLIcfg['BPFcfg_LF']['f2'] - TLIcfg['BPFcfg_LF']['f1']
    lowestFreq = f0_LF - Bw_LF / 2.0
    
    sampleT = int(round(NT * fs / lowestFreq))
    halfSampleT = sampleT // 2
    Nsamples = LFsignal.shape[0]
    
    HF_z = function_zscore_v0(HFsignal).flatten()
    LF_z = function_zscore_v0(LFsignal).flatten()
    
    # Hilbert transform for phase extraction
    LFphase = np.angle(hilbert(LF_z))
    
    s1 = LFphase[:-1]
    s2 = LFphase[1:]
    
    if LFphase_type == 'peaks':
        pulseZC = (s1 * s2 < 0) & (np.abs(s1 - s2) < np.pi)
    else:
        pulseZC = (s1 * s2 < 0) & (np.abs(s1 - s2) >= np.pi)
        
    pulseZC = np.concatenate(([False], pulseZC))
    indPeak_LF = np.where(pulseZC)[0]
    Npeak = len(indPeak_LF)
    
    if Npeak < 2:
        return {'TLI': np.nan, 'indPeak_HF': [], 'indPeak_LF': indPeak_LF}

    HFsignal_aux1 = np.abs(HF_z) if TLIcfg.get('abs_HFSignaltimeLockedHFpeaks') else HF_z
    HFsignal_aux2 = np.abs(HF_z) if TLIcfg.get('abs_HFSignaltimeLockedLFpeaks') else HF_z
    
    # Initialize integer array for indices
    indPeak_HF = np.zeros(Npeak - 1, dtype=int)
    
    for ii in range(Npeak - 1):
        segment = HFsignal_aux1[indPeak_LF[ii]:indPeak_LF[ii+1] + 1]
        indPeak_HF[ii] = np.argmax(segment) + indPeak_LF[ii]
        
    JJstart, JJend = None, None
    for jj in range(Npeak - 1):
        if (indPeak_LF[jj] - halfSampleT >= 0) and (indPeak_HF[jj] - halfSampleT >= 0):
            JJstart = jj
            break
            
    for jj in range(Npeak - 2, -1, -1):
        if (indPeak_LF[jj] + halfSampleT < Nsamples) and (indPeak_HF[jj] + halfSampleT < Nsamples):
            JJend = jj
            break
            
    if JJstart is None or JJend is None or JJstart > JJend:
        return {'TLI': np.nan, 'indPeak_HF': indPeak_HF, 'indPeak_LF': indPeak_LF}
        
    HF_locked_HF = np.zeros(2 * halfSampleT + 1)
    HF_locked_LF = np.zeros(2 * halfSampleT + 1)
    
    for jj in range(JJstart, JJend + 1):
        HF_locked_HF += HFsignal_aux1[indPeak_HF[jj] - halfSampleT : indPeak_HF[jj] + halfSampleT + 1]
        HF_locked_LF += HFsignal_aux2[indPeak_LF[jj] - halfSampleT : indPeak_LF[jj] + halfSampleT + 1]
        
    epochs_count = JJend - JJstart + 1
    HF_locked_HF /= epochs_count
    HF_locked_LF /= epochs_count
    
    denom = np.max(HF_locked_HF) - np.min(HF_locked_HF)
    numer = np.max(HF_locked_LF) - np.min(HF_locked_LF)
    
    TLI = numer / denom if denom != 0 else np.nan
    
    if TLIcfg.get('abs_HFSignaltimeLockedHFpeaks'): TLI /= 2.0
    if TLIcfg.get('abs_HFSignaltimeLockedLFpeaks'): TLI *= 2.0
    
    return {
        'TLI': TLI,
        'HFSignaltimeLockedHFpeaks': HF_locked_HF,
        'HFSignaltimeLockedLFpeaks': HF_locked_LF,
        'indPeak_HF': indPeak_HF,
        'indPeak_LF': indPeak_LF
    }


# ==========================================
# 5. Van der Pol Differential Equations
# ==========================================
def function_vanDerPolDiffEq_v0(t, z, mu, w):
    """
    Van der Pol Differential Equation.
    dz/dt = [y, mu*(1 - x^2)*y - w^2 * x]
    """
    return [z[1], mu * (1.0 - z[0]**2) * z[1] - (w**2) * z[0]]


# ==========================================
# 6. Van der Pol Solver
# ==========================================
def function_vanDerPolSolver_v0(mu_array, w, time_f, dt):
    """
    Solves the Van der Pol equation across mu parameters using Radau for stiffness.
    """
    num_points = int(round(time_f / dt)) + 1
    timeVector = np.linspace(0, time_f, num_points)
    Nmu = len(mu_array)
    signal = np.zeros((len(timeVector), Nmu))
    
    for ii, mu in enumerate(mu_array):
        # 'Radau' implicit solver handles stiff system dynamics for high mu
        sol = solve_ivp(
            fun=lambda t, z: function_vanDerPolDiffEq_v0(t, z, mu, w),
            t_span=(0, time_f),
            y0=[2.0, 1.0],
            t_eval=timeVector,
            method='Radau' if mu > 10 else 'RK45'
        )
        signal[:, ii] = sol.y[0, :]
        
    return signal


# ==========================================
# 7. Main Test Script
# ==========================================
if __name__ == "__main__":
    
    # Oscillator parameters
    f = 10.0
    w = 2.0 * np.pi * f
    mu_array = np.arange(0, 252, 2)  # 0 to 250 with step 2
    
    # Simulation settings
    dt = 0.001
    fs = 1.0 / dt
    TIME_SIM = 6.0
    SOLVER_TRANSIENT = int(fs / f)
    
    FREQ_BANDS = np.array([[1.0, 15.0], [20.0, 100.0]])
    
    # Low-Frequency BPF Configuration
    BPFcfg_LF = {
        'f1': FREQ_BANDS[0, 0],
        'f2': FREQ_BANDS[0, 1],
        'zeropadding': 0,
        'conv': 'circular',
        'causal': False,
        'fs': fs
    }
    
    # High-Frequency BPF Configuration
    BPFcfg_HF = {
        'f1': FREQ_BANDS[1, 0],
        'f2': FREQ_BANDS[1, 1],
        'zeropadding': 0,  # Fixed typo ('zeroppadding')
        'conv': 'circular',
        'causal': False,
        'fs': fs
    }
    
    # TLI Configuration
    TLIcfg = {
        'BPFcfg_LF': BPFcfg_LF,
        'BPFcfg_HF': BPFcfg_HF,
        'abs_HFSignaltimeLockedHFpeaks': False,
        'abs_HFSignaltimeLockedLFpeaks': False,
        'LFphase': 'peaks',
        'NT': 1,
        'fs': fs
    }
    
    NOISE_LEVEL = 0.1
    np.random.seed(0)  # Reproducible random noise
    
    print("Simulating Van der Pol oscillator across mu parameter space...")
    signal = function_vanDerPolSolver_v0(mu_array, w, TIME_SIM, dt)
    signal = signal[SOLVER_TRANSIENT:, :]  # Remove transient
    
    Ns = signal.shape[0]
    t = np.arange(0, Ns) / fs
    
    TLI_results = np.full(len(mu_array), np.nan)
    
    for ii in range(len(mu_array)):
        if (ii + 1) % 25 == 0 or ii == 0 or ii == len(mu_array) - 1:
            print(f"Processing: mu = {mu_array[ii]}; {ii+1} out of {len(mu_array)}.")
        
        rawSignal = signal[:, ii:ii+1]
        
        # Additive white Gaussian noise
        stdDevSignal = np.std(rawSignal, axis=0, ddof=1)
        noise = (NOISE_LEVEL * stdDevSignal) * np.random.randn(*rawSignal.shape)
        rawSignal = rawSignal + noise
        
        # Z-score normalization
        rawSignal = function_zscore_v0(rawSignal)
        
        # Edge reflection
        rawSignal_reflected = np.vstack([rawSignal[::-1], rawSignal, rawSignal[::-1]])
        
        # Filter signals
        FDFout_LF = function_FDF_v0(rawSignal_reflected, BPFcfg_LF)
        LFSIGNAL = FDFout_LF['filteredSignal']
        
        FDFout_HF = function_FDF_v0(rawSignal_reflected, BPFcfg_HF)
        HFSIGNAL = FDFout_HF['filteredSignal']
        
        # Restore length
        LFSIGNAL = LFSIGNAL[Ns:2*Ns]
        HFSIGNAL = HFSIGNAL[Ns:2*Ns]
        rawSignal_clean = rawSignal_reflected[Ns:2*Ns]
        
        # Compute TLI
        TLIout = function_TimeLockedIndex_v0(LFSIGNAL, HFSIGNAL, TLIcfg)
        TLI_results[ii] = TLIout['TLI']
        
        # Plot initial and final runs
        if ii == 0 or ii == len(mu_array) - 1:
            plt.figure(figsize=(10, 4))
            plt.plot(t, rawSignal_clean, '-k', label='Raw signal')
            plt.plot(t, LFSIGNAL, '-b', label='LF signal')
            plt.plot(t, HFSIGNAL, '-r', label='HF signal')
            plt.legend(loc='upper right')
            plt.title(f'mu / w = {mu_array[ii]/w:.3g}')
            plt.xlabel('Time [sec.]')
            plt.ylabel('Amplitude [arb. units]')
            plt.tight_layout()
            plt.show()

    # Final summary plot
    plt.figure(figsize=(8, 4.5))
    plt.plot(mu_array / w, TLI_results, ':ob', markersize=6)
    plt.grid(True)
    plt.ylim([0, 1.1])
    plt.xlabel('mu / w [arb. units]')
    plt.ylabel('TLI [arb. units]')
    plt.title('TLI as a function of parameter mu')
    plt.tight_layout()
    plt.show()