print("Script starting...")
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

def load_strain(file_path):
    with h5.File(file_path, 'r') as f:
        # strain data is stored under this path in GWOSC files
        strain = f['strain']['Strain'][:]
        # get the sampling rate and start time from metadata
        dt= f['strain']['Strain'].attrs['Xspacing']
        t0 = f['strain']['Strain'].attrs['Xstart']
    return strain, dt, t0

# load both detectors
h1_strain, dt, t0 = load_strain('data/H1_GW150914.hdf5')
l1_strain, _, _ = load_strain('data/L1_GW150914.hdf5')

# build time array
n = len(h1_strain)
t = t0 + np.arange(n) * dt

print(f"Sample rate: {1/dt:.0f} Hz")
print(f"Duration: {n * dt:.1f} seconds")
print(f"Number of samples: {n}")

# plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

ax1.plot(t - t0, h1_strain, 'b', lw=0.3)
ax1.set_ylabel('Strain')
ax1.set_title('LIGO Hanford (H1) — Raw Strain')

ax2.plot(t - t0, l1_strain, 'r', lw=0.3)
ax2.set_ylabel('Strain')
ax2.set_title('LIGO Livingston (L1) — Raw Strain')
ax2.set_xlabel('Time (seconds)')

plt.tight_layout()
plt.savefig('figures/01_raw_strain.png', dpi=150)
plt.show()