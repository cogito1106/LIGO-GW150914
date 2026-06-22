import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, welch
from scipy.signal.windows import tukey

print("Script starting...")

def load_strain(filepath):
    with h5.File(filepath, 'r') as f:
        strain = f['strain']['Strain'][:]
        dt = f['strain']['Strain'].attrs['Xspacing']
        t0 = f['strain']['Strain'].attrs['Xstart']
    return strain, dt, t0

# load data
h1_strain, dt, t0 = load_strain('data/H1_GW150914.hdf5')
l1_strain, _, _  = load_strain('data/L1_GW150914.hdf5')
fs = 1.0 / dt

print(f"Sample rate: {fs:.0f} Hz")

# compute PSD using Welch's method
freqs_h1, psd_h1 = welch(h1_strain, fs=fs, nperseg=4*fs)
freqs_l1, psd_l1 = welch(l1_strain, fs=fs, nperseg=4*fs)

# plot the PSD
fig, ax = plt.subplots(figsize=(10, 5))
ax.loglog(freqs_h1, np.sqrt(psd_h1), 'b', label='H1', alpha=0.8)
ax.loglog(freqs_l1, np.sqrt(psd_l1), 'r', label='L1', alpha=0.8)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Strain noise [1/√Hz]')
ax.set_title('LIGO Noise Amplitude Spectral Density')
ax.legend()
ax.set_xlim([10, 2000])
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/02_psd.png', dpi=150)
plt.show()

def whiten(strain, psd, freqs, dt):
    n = len(strain)

    # taper edges to zero to prevent edge artefacts
    window = tukey(n, alpha=0.1)
    strain_windowed = strain * window

    # FFT
    strain_fft = np.fft.rfft(strain_windowed)
    fft_freqs = np.fft.rfftfreq(n, d=dt)

    # interpolate PSD onto FFT frequency grid
    psd_interp = np.interp(fft_freqs, freqs, psd)

    # avoid division by zero at DC
    psd_interp[0] = psd_interp[1]

    # divide by amplitude spectral density
    whitened_fft = strain_fft / np.sqrt(psd_interp)

    # inverse FFT back to time domain
    whitened = np.fft.irfft(whitened_fft, n=n)

    return whitened

def bandpass(strain, lowcut, highcut, fs, order=4):
    nyq = fs / 2
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, strain)

# whiten both detectors
h1_white = whiten(h1_strain, psd_h1, freqs_h1, dt)
l1_white = whiten(l1_strain, psd_l1, freqs_l1, dt)

# bandpass filter 35-350 Hz
h1_bp = bandpass(h1_white, 35, 350, fs)
l1_bp = bandpass(l1_white, 35, 350, fs)

# build time array
n = len(h1_strain)
t = t0 + np.arange(n) * dt

# plot whitened and bandpassed strain
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

ax1.plot(t - t0, h1_bp, 'b', lw=0.5)
ax1.set_ylabel('Whitened Strain')
ax1.set_title('LIGO Hanford (H1) — Whitened & Bandpassed')
ax1.set_xlim([0, 32])

ax2.plot(t - t0, l1_bp, 'r', lw=0.5)
ax2.set_ylabel('Whitened Strain')
ax2.set_title('LIGO Livingston (L1) — Whitened & Bandpassed')
ax2.set_xlabel('Time (seconds)')
ax2.set_xlim([0, 32])

plt.tight_layout()
plt.savefig('figures/03_whitened.png', dpi=150)
plt.show()

print("Whitening done.")