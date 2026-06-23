# LIGO GW150914 — Gravitational Wave Detection from Scratch

Reproducing the first-ever gravitational wave detection (GW150914) using matched filtering on real LIGO strain data from the Gravitational Wave Open Science Center (GWOSC).

Built from scratch in Python — no black-box GW libraries used.

---

## What this project does

On 14 September 2015, LIGO detected gravitational waves from two merging black holes (36 M☉ and 29 M☉) located ~1.3 billion light years away. This project reproduces that detection step by step:

1. **Raw strain loading** — download and inspect real LIGO H1/L1 strain data
2. **Noise characterisation** — estimate the power spectral density via Welch's method
3. **Whitening & bandpassing** — flatten the noise floor across 35–350 Hz
4. **Template generation** — compute a post-Newtonian inspiral waveform from GR
5. **Matched filtering** — cross-correlate template with data in the frequency domain
6. **Detection** — recover SNR ~29 spike at the correct merger time in both detectors
7. **Visualisation** — spectrogram, zoomed strain, and template overlay

---

## Results

| Quantity | This project | LIGO published |
|---|---|---|
| Peak SNR (H1) | ~29 | ~24 |
| Inter-detector time delay | 7.1 ms | 6.9 ms |
| Merger time (in segment) | ~15.4 s | 15.4 s |

---

## Pipeline
---

## Key results

![SNR](figures/05_snr.png)
![Spectrogram](figures/07_spectrogram.png)
![Template Overlay](figures/08_template_overlay.png)

---

## Setup

```bash
git clone https://github.com/cogito1106/LIGO-GW150914.git
cd LIGO-GW150914
pip install -r requirements.txt
python src/download_data.py   # download LIGO data (~400MB)
```

Then run the notebooks in order.

---

## Physics background

Gravitational waves are ripples in spacetime produced by accelerating masses. For a compact binary system, the strain amplitude is:

h ~ (G Mc / c²)^(5/3) (π f)^(2/3) / (c² d)

where Mc is the chirp mass and d is the luminosity distance. The matched filter cross-correlates a theoretical template h(t) with the detector strain s(t), weighted by the noise PSD S(f):

SNR(t) = |∫ s̃(f) h̃*(f) / S(f) e^(2πift) df| / σ

The SNR peaks when the template best matches the signal buried in the data.

---

## References

- Abbott et al. (2016) — [GW150914 detection paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.061102)
- [GWOSC](https://gwosc.org) — Gravitational Wave Open Science Center
- Maggiore (2007) — Gravitational Waves: Theory and Experiments