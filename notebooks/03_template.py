import numpy as np
import matplotlib.pyplot as plt

print("Script starting...")

# GW150914 best-fit parameters from the LIGO paper
M1 = 36.0  # solar masses, heavier black hole
M2 = 29.0  # solar masses, lighter black hole

# constants
G = 6.674e-11      # gravitational constant
c = 3e8            # speed of light
M_sun = 1.989e30   # solar mass in kg

def compute_chirp_mass(m1, m2):
    """Chirp mass — the combination of masses that GWs are most sensitive to"""
    m1_kg = m1 * M_sun
    m2_kg = m2 * M_sun
    Mc = (m1_kg * m2_kg)**(3/5) / (m1_kg + m2_kg)**(1/5)
    return Mc

def generate_inspiral_template(Mc, f_low, f_high, dt):
    # time to coalescence from starting frequency
    def t_coal(f, Mc):
        return (5/256) * (c**3 / (G * Mc))**(5/3) * (np.pi * f)**(-8/3)

    t_start = t_coal(f_low, Mc)

    # time array up to coalescence
    t = np.arange(-t_start, -dt, dt)
    tau = -t  # time remaining until coalescence

    # GW frequency: f = (1/pi) * (5/(256*tau))^(3/8) * (c^3/(G*Mc))^(5/8)
    f_inst = (1.0 / np.pi) * (5.0 / (256.0 * tau))**(3.0/8.0) * (c**3 / (G * Mc))**(5.0/8.0)

    # keep only samples below f_high
    mask = f_inst < f_high
    t = t[mask]
    tau = tau[mask]
    f_inst = f_inst[mask]

    # phase by integrating frequency
    phase = 2 * np.pi * np.cumsum(f_inst) * dt

    # amplitude
    amplitude = tau**(-1/4)
    amplitude /= amplitude.max()

    h = amplitude * np.cos(phase)

    return t, h, f_inst

# compute chirp mass
Mc = compute_chirp_mass(M1, M2)
print(f"Chirp mass: {Mc/M_sun:.2f} solar masses")

# generate template
dt = 1/4096
t, h, f_inst = generate_inspiral_template(Mc, f_low=20.0, f_high=500.0, dt=dt)

print(f"Template duration: {len(t)*dt:.2f} seconds")

# plot waveform
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

ax1.plot(t, h, 'k', lw=0.5)
ax1.set_xlabel('Time to coalescence (s)')
ax1.set_ylabel('Strain (arbitrary units)')
ax1.set_title('Post-Newtonian Inspiral Waveform — GW150914 Parameters')
ax1.grid(True, alpha=0.3)

ax2.plot(t, f_inst, 'purple', lw=1)
ax2.set_xlabel('Time to coalescence (s)')
ax2.set_ylabel('Frequency (Hz)')
ax2.set_title('Instantaneous Frequency — The Chirp')
ax2.set_ylim([0, 500])
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/04_template.png', dpi=150)
plt.show()

print("Template generated.")