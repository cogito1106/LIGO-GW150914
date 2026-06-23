import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, welch, spectrogram
from scipy.signal.windows import tukey

print("Script starting...")

# constants
G = 6.674e-11
c = 3e8
M_sun = 1.989e30
M1, M2 = 36.0, 29.0

def load_strain(filepath):
    with h5.File(filepath, 'r') as f:
        strain = f['strain']['Strain'][:]
        dt = f['strain']['Strain'].attrs['Xspacing']
        t0 = f['strain']['Strain'].attrs['Xstart']
    return strain, dt, t0

def whiten(strain, psd, freqs, dt):
    n = len(strain)
    window = tukey(n, alpha=0.1)
    strain_fft = np.fft.rfft(strain * window)
    fft_freqs = np.fft.rfftfreq(n, d=dt)
    psd_interp = np.interp(fft_freqs, freqs, psd)
    psd_interp[0] = psd_interp[1]
    return np.fft.irfft(strain_fft / np.sqrt(psd_interp), n=n)

def bandpass(strain, lowcut, highcut, fs, order=4):
    nyq = fs / 2
    b, a = butter(order, [lowcut/nyq, highcut/nyq], btype='band')
    return filtfilt(b, a, strain)

def compute_chirp_mass(m1, m2):
    m1_kg, m2_kg = m1 * M_sun, m2 * M_sun
    return (m1_kg * m2_kg)**(3/5) / (m1_kg + m2_kg)**(1/5)

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
    return t, amplitude * np.cos(phase)

# load data
h1_strain, dt, t0 = load_strain('data/H1_GW150914.hdf5')
fs = 1.0 / dt
n = len(h1_strain)
t = t0 + np.arange(n) * dt
t_rel = t - t0  # time relative to start

# PSD and whitening
freqs_h1, psd_h1 = welch(h1_strain, fs=fs, nperseg=4*fs)
h1_white = bandpass(whiten(h1_strain, psd_h1, freqs_h1, dt), 35, 350, fs)

# known merger time in this 32s segment (GPS 1126259462.4 → ~15.4s into segment)
t_merger = 15.4

# ---- PLOT 1: zoomed whitened strain around merger ----
zoom_window = 0.4  # seconds either side
mask_zoom = np.abs(t_rel - t_merger) < zoom_window

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(t_rel[mask_zoom] - t_merger, h1_white[mask_zoom], 'b', lw=1)
ax.set_xlabel('Time relative to merger (s)')
ax.set_ylabel('Whitened Strain')
ax.set_title('GW150914 — Zoomed Whitened Strain around Merger (H1)')
ax.axvline(x=0, color='r', linestyle='--', alpha=0.5, label='Merger time')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/06_zoomed_strain.png', dpi=150)
plt.show()

# ---- PLOT 2: spectrogram ----
# zoom to 4 seconds around merger
mask_spec = np.abs(t_rel - t_merger) < 2.0
h1_zoom = h1_strain[mask_spec]

f_spec, t_spec, Sxx = spectrogram(
    h1_zoom, fs=fs,
    nperseg=256,
    noverlap=240,
    window='hann'
)

fig, ax = plt.subplots(figsize=(10, 5))
ax.pcolormesh(
    t_spec - 2.0, f_spec, np.log10(Sxx),
    shading='gouraud', cmap='inferno',
    vmin=-46, vmax=-40
)
ax.set_ylim([20, 400])
ax.set_xlabel('Time relative to merger (s)')
ax.set_ylabel('Frequency (Hz)')
ax.set_title('GW150914 Spectrogram — H1 (The Chirp)')
ax.axvline(x=0, color='w', linestyle='--', alpha=0.7, label='Merger time')
ax.legend()
plt.tight_layout()
plt.savefig('figures/07_spectrogram.png', dpi=150)
plt.show()

# ---- PLOT 3: template overlay ----
Mc = compute_chirp_mass(M1, M2)
t_tmpl, h_tmpl = generate_template(Mc, f_low=20.0, f_high=500.0, dt=dt)

# normalise template to match whitened strain amplitude
h_tmpl_norm = h_tmpl * (np.max(np.abs(h1_white[mask_zoom])) / np.max(np.abs(h_tmpl))) * 0.6

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(t_rel[mask_zoom] - t_merger, h1_white[mask_zoom], 'b', lw=1, label='H1 whitened strain', alpha=0.8)
ax.plot(t_tmpl, h_tmpl_norm, 'r', lw=1.5, label='PN template', alpha=0.9)
ax.set_xlabel('Time relative to merger (s)')
ax.set_ylabel('Whitened Strain')
ax.set_title('GW150914 — Whitened Strain vs Template Overlay (H1)')
ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/08_template_overlay.png', dpi=150)
plt.show()

print("All visualisations saved.")