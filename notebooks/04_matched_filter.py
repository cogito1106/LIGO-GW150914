import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, welch
from scipy.signal.windows import tukey

print("Script starting...")

# constants
G = 6.674e-11
c = 3e8
M_sun = 1.989e30

# GW150914 parameters
M1 = 36.0
M2 = 29.0

def load_strain(filepath):
    with h5.File(filepath, 'r') as f:
        strain = f['strain']['Strain'][:]
        dt = f['strain']['Strain'].attrs['Xspacing']
        t0 = f['strain']['Strain'].attrs['Xstart']
    return strain, dt, t0

def compute_chirp_mass(m1, m2):
    m1_kg = m1 * M_sun
    m2_kg = m2 * M_sun
    return (m1_kg * m2_kg)**(3/5) / (m1_kg + m2_kg)**(1/5)

def whiten(strain, psd, freqs, dt):
    n = len(strain)
    window = tukey(n, alpha=0.1)
    strain_fft = np.fft.rfft(strain * window)
    fft_freqs = np.fft.rfftfreq(n, d=dt)
    psd_interp = np.interp(fft_freqs, freqs, psd)
    psd_interp[0] = psd_interp[1]
    whitened_fft = strain_fft / np.sqrt(psd_interp)
    return np.fft.irfft(whitened_fft, n=n)

def bandpass(strain, lowcut, highcut, fs, order=4):
    nyq = fs / 2
    b, a = butter(order, [lowcut/nyq, highcut/nyq], btype='band')
    return filtfilt(b, a, strain)

def generate_template(Mc, f_low, f_high, dt):
    def t_coal(f):
        return (5/256) * (c**3 / (G * Mc))**(5/3) * (np.pi * f)**(-8/3)
    t = np.arange(-t_coal(f_low), -dt, dt)
    tau = -t
    f_inst = (1.0/np.pi) * (5.0/(256.0*tau))**(3/8) * (c**3/(G*Mc))**(5/8)
    mask = f_inst < f_high
    t, tau, f_inst = t[mask], tau[mask], f_inst[mask]
    phase = 2 * np.pi * np.cumsum(f_inst) * dt
    amplitude = tau**(-1/4)
    amplitude /= amplitude.max()
    return amplitude * np.cos(phase)

def matched_filter(data, template, dt, psd, freqs):
    n = len(data)
    
    # zero pad template to same length as data
    template_padded = np.zeros(n)
    template_padded[:len(template)] = template
    
    # FFT both
    data_fft = np.fft.rfft(data)
    template_fft = np.fft.rfft(template_padded)
    fft_freqs = np.fft.rfftfreq(n, d=dt)
    
    # interpolate PSD onto FFT grid
    psd_interp = np.interp(fft_freqs, freqs, psd)
    psd_interp[0] = psd_interp[1]
    
    # matched filter in frequency domain:
    # cross correlate data with template, weighted by noise PSD
    mf_fft = (data_fft * np.conj(template_fft)) / psd_interp
    
    # inverse FFT to get SNR time series
    snr_complex = np.fft.irfft(mf_fft, n=n)
    
    # normalise to get SNR
    sigma = np.sqrt(2 * np.sum(np.abs(template_fft)**2 / psd_interp) * dt / n**2)
    snr = np.abs(snr_complex) / sigma
    
    return snr

# load data
h1_strain, dt, t0 = load_strain('data/H1_GW150914.hdf5')
l1_strain, _, _  = load_strain('data/L1_GW150914.hdf5')
fs = 1.0 / dt
n = len(h1_strain)
t = t0 + np.arange(n) * dt

print(f"Sample rate: {fs:.0f} Hz, Duration: {n*dt:.1f}s")

# compute PSDs
freqs_h1, psd_h1 = welch(h1_strain, fs=fs, nperseg=4*fs)
freqs_l1, psd_l1 = welch(l1_strain, fs=fs, nperseg=4*fs)

# whiten and bandpass
h1_white = bandpass(whiten(h1_strain, psd_h1, freqs_h1, dt), 35, 350, fs)
l1_white = bandpass(whiten(l1_strain, psd_l1, freqs_l1, dt), 35, 350, fs)

# generate template
Mc = compute_chirp_mass(M1, M2)
template = generate_template(Mc, f_low=20.0, f_high=500.0, dt=dt)
print(f"Template length: {len(template)*dt:.2f}s")

# run matched filter on both detectors
snr_h1 = matched_filter(h1_white, template, dt, psd_h1, freqs_h1)
snr_l1 = matched_filter(l1_white, template, dt, psd_l1, freqs_l1)

# plot SNR time series
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

ax1.plot(t - t0, snr_h1, 'b', lw=0.8)
ax1.axhline(y=8, color='k', linestyle='--', alpha=0.5, label='SNR=8 threshold')
ax1.set_ylabel('SNR')
ax1.set_title('Matched Filter SNR — LIGO Hanford (H1)')
ax1.legend()
ax1.set_xlim([0, 32])

ax2.plot(t - t0, snr_l1, 'r', lw=0.8)
ax2.axhline(y=8, color='k', linestyle='--', alpha=0.5, label='SNR=8 threshold')
ax2.set_ylabel('SNR')
ax2.set_title('Matched Filter SNR — LIGO Livingston (L1)')
ax2.set_xlabel('Time (seconds)')
ax2.legend()
ax2.set_xlim([0, 32])

plt.tight_layout()
plt.savefig('figures/05_snr.png', dpi=150)
plt.show()

# find peak SNR times
peak_h1 = t[np.argmax(snr_h1)] - t0
peak_l1 = t[np.argmax(snr_l1)] - t0
print(f"H1 peak SNR: {snr_h1.max():.1f} at t={peak_h1:.3f}s")
print(f"L1 peak SNR: {snr_l1.max():.1f} at t={peak_l1:.3f}s")
print(f"Time delay between detectors: {(peak_h1 - peak_l1)*1000:.1f} ms")