# ConformalAntennaSimulation

FDTD electromagnetic simulation (openEMS) of a 1.575 GHz microstrip patch antenna, evaluated flat and conformally bent (r=90mm), with S11, impedance, directivity, and radiation pattern analysis.

---

## Overview

This project designs and simulates a microstrip patch antenna targeting the **L-band (1.575 GHz)**, evaluated in two configurations:

- **Flat** — standard patch antenna on a flat FR4 substrate (baseline)
- **Bent / conformal** — the same design curved around a 90mm-radius cylindrical surface, representing mounting on a curved surface (e.g. a helmet)

The core objective: quantify how bending a patch antenna onto a curved surface affects its resonant frequency, impedance matching, and radiation pattern — a practical question for any wearable or conformal antenna application.

Simulations were performed using **openEMS**, a free, open-source electromagnetic solver that numerically solves Maxwell's equations via the **FDTD (Finite-Difference Time-Domain)** method.

---

## Design Parameters

| Parameter | Value |
|---|---|
| Target frequency | 1.575 GHz (L-band) |
| Substrate | FR4, εr = 4.4, thickness 1.6mm |
| Patch width | 44.0 mm |
| Patch length | 34.2 mm |
| Feed type | Probe-fed |
| Feed offset | -5.2 mm from patch center |
| Feed impedance | 50 Ω |
| Bend radius (conformal case) | 90 mm |

---

## Results

| Metric | Flat | Bent (r = 90mm) |
|---|---|---|
| Resonant frequency | 1.585 GHz | 1.555 GHz |
| S11 at resonance | -6.77 dB | -6.43 dB |
| Directivity | — | 5.8 dBi |
| Radiation efficiency | — | 28.5% |
| Radiation pattern | Broad forward lobe | Broad forward lobe, reshaped by curvature |

**Key finding:** bending the patch onto a 90mm-radius curve shifted resonance by **~30 MHz (1.9%)**. This confirms that curvature measurably changes a patch antenna's effective electrical length and must be accounted for at design time — a flat-tuned design cannot simply be mounted onto a curved surface without retuning.

---

## Repository Contents

| File | Description |
|---|---|
| `L_band_patch_antenna_1575MHz.py` | Flat patch antenna simulation |
| `Bent_L_band_patch_antenna_helmet.py` | Bent/conformal patch antenna simulation |
| `Simulation_Summary.md` | Full write-up: methodology, results, and known limitations |

---

## Requirements

- Python 3.x
- [openEMS](https://docs.openems.de/) + CSXCAD (Python bindings)
- `numpy`, `matplotlib`, `h5py`

On Windows, each script requires a DLL loading fix before importing openEMS/CSXCAD:

```python
import os
os.add_dll_directory(r"path/to/your/openEMS/install")
```

---

## Status / Roadmap

- [x] Flat patch design — simulated and tuned
- [x] Bent/conformal patch design — simulated and compared against flat baseline
- [ ] Second band (UHF-class, ~900 MHz) — dimensions calculated (101.4mm × 79.3mm on FR4), not yet simulated
- [ ] Improve impedance matching beyond -10dB (currently ~-6.5dB on both configurations)
- [ ] Full doubly-curved (dome) surface model, beyond the current cylindrical approximation
- [ ] Time-domain near-field propagation animation

---

## Acknowledgements

Built on the official openEMS tutorial scripts (`Simple_Patch_Antenna` and `Bent_Patch_Antenna`) from the [openEMS-Project](https://github.com/thliebig/openEMS-Project).