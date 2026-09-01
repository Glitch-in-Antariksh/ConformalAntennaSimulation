# -*- coding: utf-8 -*-
"""
Helmet Conformal Antenna Project - SIH
L-band Patch Antenna Simulation (target: 1575 MHz)

Adapted from the official openEMS "Simple Patch Antenna" tutorial:
https://raw.githubusercontent.com/thliebig/openEMS/master/python/Tutorials/Simple_Patch_Antenna.py

Substrate: FR4 (epsR=4.4), thickness 1.6mm
Patch dimensions calculated using standard transmission-line model formulas
for f0 = 1575 MHz (L-band, matches GPS L1 - used as our video/camera link band)
"""

### Import Libraries -- NOTE: DLL fix required on Windows before these imports work
import os
os.add_dll_directory(r"D:\ANTARIKSH\openEMS\openEMS")  # <-- adjust if your path differs

import tempfile
import numpy as np
import matplotlib.pyplot as plt

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Lband_Patch')

post_proc_only = False

# --- Patch dimensions (RETUNED again after 2nd sim run) ---
# 2nd attempt (42.3 x 32.9mm) resonated at ~1.65 GHz, overshooting the
# 1.575 GHz target. Interpolating between both data points (45.13mm->1.15GHz,
# 32.9mm->1.65GHz) to land closer to target this time.
patch_width  = 44.0   # mm, x-direction
patch_length = 34.2   # mm, y-direction (resonant length)

# --- Substrate setup (FR4) ---
substrate_epsR   = 4.4
tan_delta        = 0.02   # typical FR4 loss tangent
f0 = 1.575e9   # center frequency: L-band, matches your GPS-style target
fc = 0.5e9     # 20dB corner frequency (keeps sim focused near L-band)

substrate_kappa  = tan_delta * 2*np.pi*f0 * EPS0*substrate_epsR
substrate_width  = 90    # mm, some margin around the patch
substrate_length = 80    # mm
substrate_thickness = 1.6   # mm, standard FR4 thickness
substrate_cells = 4

# --- Feeding position ---
# Inset feed position (x-direction), scaled proportionally from the original
# tutorial's -6mm/40mm ratio. This is a STARTING GUESS -- if S11 shows poor
# matching (not a deep enough dip near 1575 MHz), try shifting this value
# in a few mm steps and re-running.
feed_pos = -5.2   # mm (scaled proportionally with new patch size)
feed_R = 50        # feed resistance (ohms)

# size of the simulation "air box" surrounding the structure
SimBox = np.array([220, 220, 160])

### FDTD setup
FDTD = openEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite( f0, fc )
FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)
mesh_res = C0/(f0+fc)/1e-3/20

### Generate properties, primitives and mesh-grid
mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2])
mesh.AddLine('z', [-SimBox[2]/3, SimBox[2]*2/3])

# create patch
patch = CSX.AddMetal('patch')
start = [-patch_width/2, -patch_length/2, substrate_thickness]
stop  = [ patch_width/2 , patch_length/2, substrate_thickness]
patch.AddBox(priority=10, start=start, stop=stop)
FDTD.AddEdges2Grid(dirs='xy', properties=patch, metal_edge_res=mesh_res/2)

# create substrate
substrate = CSX.AddMaterial('substrate', epsilon=substrate_epsR, kappa=substrate_kappa)
start = [-substrate_width/2, -substrate_length/2, 0]
stop  = [ substrate_width/2,  substrate_length/2, substrate_thickness]
substrate.AddBox(priority=0, start=start, stop=stop)

mesh.AddLine('z', np.linspace(0, substrate_thickness, substrate_cells+1))

# create ground plane
gnd = CSX.AddMetal('gnd')
start[2] = 0
stop[2]  = 0
gnd.AddBox(start, stop, priority=10)
FDTD.AddEdges2Grid(dirs='xy', properties=gnd)

# feed / port
start = [feed_pos, 0, 0]
stop  = [feed_pos, 0, substrate_thickness]
port = FDTD.AddLumpedPort(1, feed_R, start, stop, 'z', 1.0, priority=5, edges2grid='xy')

mesh.SmoothMeshLines('all', mesh_res, 1.4)

# nf2ff recording box (for far-field radiation pattern)
nf2ff = FDTD.CreateNF2FFBox()

### View the 3D geometry (opens a separate window; close it to continue)
# Set to False to skip the viewer and go straight to simulation.
view_geometry = True
if view_geometry:
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX_file = os.path.join(Sim_Path, 'lband_patch.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

### Run the simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-processing and plotting
f = np.linspace(max(1e9, f0-fc), f0+fc, 401)
port.CalcPort(Sim_Path, f)
s11 = port.uf_ref/port.uf_inc
s11_dB = 20.0*np.log10(np.abs(s11))

fig, axis = plt.subplots(num="S11", tight_layout=True)
axis.plot(f/1e9, s11_dB, 'k-', linewidth=2, label='S11')
axis.axvline(f0/1e9, color='b', linestyle=':', label='Target 1.575 GHz')
axis.grid()
axis.set_xmargin(0)
axis.set_xlabel('Frequency (GHz)')
axis.set_ylabel('S-Parameter (dB)')
axis.set_title("Input Matching - L-band Patch")
axis.legend()

idx = np.argmin(s11_dB)   # just take the best match, whatever it is
f_res = f[idx]
best_match_dB = s11_dB[idx]
print(f'Best match found at {f_res/1e9:.4f} GHz, S11 = {best_match_dB:.2f} dB')
if best_match_dB > -10:
    print('Note: match is weaker than -10dB (not ideal), but proceeding '
          'with far-field calc anyway. Consider further dimension tuning.')
if True:
    theta = np.arange(-180.0, 180.0, 2.0)
    phi   = [0., 90.]
    nf2ff_res = nf2ff.CalcNF2FF(Sim_Path, f_res, theta, phi, center=[0,0,1e-3])

    E_norm = 20.0*np.log10(nf2ff_res.E_norm[0]/np.max(nf2ff_res.E_norm[0])) + 10.0*np.log10(nf2ff_res.Dmax[0])
    fig, axis = plt.subplots(num="Pattern", tight_layout=True)
    axis.plot(theta, np.squeeze(E_norm[:,0]), 'k-', linewidth=2, label='xz-plane')
    axis.plot(theta, np.squeeze(E_norm[:,1]), 'r--', linewidth=2, label='yz-plane')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Theta (deg)')
    axis.set_ylabel('Directivity (dBi)')
    axis.set_title(f'Radiation Pattern @ {f_res/1e9:.3f} GHz')
    axis.legend()

# --- Full 3D radiation pattern ---
# Sweep the full sphere (not just 2 flat slices) and plot as a 3D lobe shape.
if True:
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (enables 3D projection)

    theta_3d = np.arange(0.0, 180.0, 4.0)      # coarser steps = faster, still smooth
    phi_3d   = np.arange(-180.0, 180.0, 4.0)
    nf2ff_3d = nf2ff.CalcNF2FF(Sim_Path, f_res, theta_3d, phi_3d, center=[0,0,1e-3])

    E_far = nf2ff_3d.E_norm[0] / np.max(nf2ff_3d.E_norm[0])  # normalize 0-1

    THETA, PHI = np.meshgrid(theta_3d*np.pi/180, phi_3d*np.pi/180, indexing='ij')

    # convert far-field magnitude + angles into 3D Cartesian "lobe" shape
    X = E_far * np.sin(THETA) * np.cos(PHI)
    Y = E_far * np.sin(THETA) * np.sin(PHI)
    Z = E_far * np.cos(THETA)

    fig = plt.figure(num="3D Radiation Pattern", figsize=(8, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, facecolors=plt.cm.jet(E_far),
                            rstride=1, cstride=1, linewidth=0, antialiased=True)
    ax.set_title(f'3D Radiation Pattern @ {f_res/1e9:.3f} GHz')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

Zin = port.uf_tot/port.if_tot
fig, axis = plt.subplots(num="Zin", tight_layout=True)
axis.plot(f/1e9, np.real(Zin), 'k-', linewidth=2, label='Re{Zin}')
axis.plot(f/1e9, np.imag(Zin), 'r--', linewidth=2, label='Im{Zin}')
axis.grid()
axis.set_xmargin(0)
axis.set_xlabel('Frequency (GHz)')
axis.set_ylabel('Zin (Ohm)')
axis.set_title("Input Impedance")
axis.legend()

plt.show()