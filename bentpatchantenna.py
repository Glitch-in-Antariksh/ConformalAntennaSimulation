# -*- coding: utf-8 -*-
"""
Helmet Conformal Antenna Project - SIH
BENT / CONFORMAL L-band Patch Antenna Simulation (target: 1575 MHz)

Adapted from the official openEMS "Bent Patch Antenna" tutorial:
https://docs.openems.de/python/openEMS/Tutorials/Bent_Patch_Antenna.html

This is the CONFORMAL version of L_band_patch_antenna_1575MHz.py -- same
patch design, but wrapped around a curved (cylindrical) surface instead of
sitting flat. patch_radius represents the curvature of the mounting surface
(here: an approximate helmet curvature radius).

Compare this script's results against the flat version to show how bending
affects resonance and radiation pattern -- this comparison is the core
proof-of-concept for the "conformal helmet antenna" problem statement.
"""

### Import Libraries -- DLL fix required on Windows before these imports work
import os
os.add_dll_directory(r"D:\ANTARIKSH\openEMS\openEMS")  # <-- adjust if your path differs

import tempfile
import numpy as np
import matplotlib.pyplot as plt

from CSXCAD import CSXCAD
from openEMS.openEMS import openEMS
from openEMS.physical_constants import *

### General parameter setup
Sim_Path = os.path.join(tempfile.gettempdir(), 'Bent_Lband_Patch')

post_proc_only = False

unit = 1e-3  # all lengths in mm

f0 = 1.575e9   # center frequency: L-band target (same as flat version)
lambda0 = round(C0/f0/unit)  # wavelength in mm, just for reference
fc = 0.5e9     # 20dB corner frequency

# --- Patch dimensions (same tuned values as the flat version) ---
patch_width  = 44.0   # mm, resonant width (alpha/angular direction when bent)
patch_length = 34.2   # mm, patch length (z-direction, along the bend axis)

# --- Curvature radius: represents the helmet mounting surface ---
# Adult head circumference is roughly 550-580mm -> radius ~90mm.
# This is the KEY new parameter vs. the flat version.
patch_radius = 90     # mm

# --- Substrate setup (FR4, same as flat version) ---
substrate_epsR   = 4.4
tan_delta        = 0.02
substrate_kappa  = tan_delta * 2*np.pi*f0 * EPS0*substrate_epsR
substrate_width  = 74     # mm, some margin around the patch
substrate_length = 64     # mm
substrate_thickness = 1.6  # mm
substrate_cells = 4

# --- Feeding ---
feed_pos = -5.2   # mm (same tuned value as flat version)
feed_R = 50

# size of the simulation domain (cylindrical: radial extent + height)
SimBox_rad    = 2*100
SimBox_height = 1.5*200

### FDTD setup -- CoordSystem=1 activates CYLINDRICAL coordinates
FDTD = openEMS(CoordSystem=1, EndCriteria=1e-4)
FDTD.SetGaussExcite(f0, fc)
FDTD.SetBoundaryCond(['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'])

# init a cylindrical mesh
CSX = CSXCAD.ContinuousStructure(CoordSystem=1)
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

### Geometry in cylindrical coordinates: [radius, angle(rad), z]
# convert linear widths into angular widths at the given radius
patch_ang_width  = patch_width/(patch_radius+substrate_thickness)
substr_ang_width = substrate_width/patch_radius
feed_angle       = feed_pos/patch_radius

# create patch
patch = CSX.AddMetal('patch')
start = [patch_radius+substrate_thickness, -patch_ang_width/2, -patch_length/2]
stop  = [patch_radius+substrate_thickness,  patch_ang_width/2,  patch_length/2]
patch.AddBox(priority=10, start=start, stop=stop)
FDTD.AddEdges2Grid(dirs='all', properties=patch)

# create substrate
substrate = CSX.AddMaterial('substrate', epsilon=substrate_epsR, kappa=substrate_kappa)
start = [patch_radius,                     -substr_ang_width/2, -substrate_length/2]
stop  = [patch_radius+substrate_thickness,  substr_ang_width/2,  substrate_length/2]
substrate.AddBox(start=start, stop=stop)
FDTD.AddEdges2Grid(dirs='all', properties=substrate)

# save current density on the patch (bonus: can visualize current flow later)
jt_patch = CSX.AddDump('Jt_patch', dump_type=3, file_type=1)
start = [patch_radius+substrate_thickness, -substr_ang_width/2, -substrate_length/2]
stop  = [patch_radius+substrate_thickness, +substr_ang_width/2,  substrate_length/2]
jt_patch.AddBox(start=start, stop=stop)

# create ground plane (curved, same radius as substrate's inner surface)
gnd = CSX.AddMetal('gnd')
start = [patch_radius, -substr_ang_width/2, -substrate_length/2]
stop  = [patch_radius, +substr_ang_width/2, +substrate_length/2]
gnd.AddBox(priority=10, start=start, stop=stop)
FDTD.AddEdges2Grid(dirs='all', properties=gnd)

# feed / port (radial direction 'r', since energy feeds through the substrate thickness)
start = [patch_radius,                     feed_angle, 0]
stop  = [patch_radius+substrate_thickness, feed_angle, 0]
port = FDTD.AddLumpedPort(1, feed_R, start, stop, 'r', 1.0, priority=50, edges2grid='all')

### Finalize the mesh
mesh.AddLine('r', patch_radius+np.array([-20, SimBox_rad]))
mesh.AddLine('a', [-0.75*np.pi, 0.75*np.pi])
mesh.AddLine('z', [-SimBox_height/2, SimBox_height/2])
mesh.AddLine('r', patch_radius+np.linspace(0, substrate_thickness, substrate_cells))

max_res = C0/(f0+fc)/unit/20
max_ang = max_res/(SimBox_rad+patch_radius)
mesh.SmoothMeshLines(0, max_res, 1.4)
mesh.SmoothMeshLines(1, max_ang, 1.4)
mesh.SmoothMeshLines(2, max_res, 1.4)

nf2ff = FDTD.CreateNF2FFBox()

### View the 3D geometry (opens a separate window; close it to continue)
view_geometry = True
if view_geometry:
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX_file = os.path.join(Sim_Path, 'bent_lband_patch.xml')
    CSX.Write2XML(CSX_file)
    from CSXCAD import AppCSXCAD_BIN
    os.system(AppCSXCAD_BIN + ' "{}"'.format(CSX_file))

### Run the simulation
if not post_proc_only:
    FDTD.Run(Sim_Path, cleanup=True)

### Post-processing
f = np.linspace(max(1e9, f0-fc), f0+fc, 401)
port.CalcPort(Sim_Path, f)
Zin = port.uf_tot / port.if_tot
s11 = port.uf_ref/port.uf_inc
s11_dB = 20.0*np.log10(np.abs(s11))

fig, axis = plt.subplots(num="S11 (bent)", tight_layout=True)
axis.plot(f/1e9, s11_dB, 'k-', linewidth=2, label='S11 (bent)')
axis.axvline(f0/1e9, color='b', linestyle=':', label='Target 1.575 GHz')
axis.grid()
axis.set_xmargin(0)
axis.set_xlabel('Frequency (GHz)')
axis.set_ylabel('S-Parameter (dB)')
axis.set_title(f"Input Matching - Bent Patch (radius={patch_radius}mm)")
axis.legend()

P_in = 0.5*np.real(port.uf_tot * np.conj(port.if_tot))  # feed power

fig, axis = plt.subplots(num="Zin (bent)", tight_layout=True)
axis.plot(f/1e6, np.real(Zin), 'k-', linewidth=2, label='Re{Zin}')
axis.plot(f/1e6, np.imag(Zin), 'r--', linewidth=2, label='Im{Zin}')
axis.grid()
axis.set_xmargin(0)
axis.set_xlabel('Frequency (MHz)')
axis.set_ylabel('Impedance (Ohm)')
axis.set_title('Feed Point Impedance (bent)')
axis.legend()

# Always compute far-field at best match (same approach as flat version)
idx = np.argmin(s11_dB)
f_res = f[idx]
best_match_dB = s11_dB[idx]
print(f'[BENT] Best match found at {f_res/1e9:.4f} GHz, S11 = {best_match_dB:.2f} dB')
if best_match_dB > -10:
    print('Note: match weaker than -10dB, proceeding with far-field calc anyway.')

theta = np.arange(-180.0, 180.0, 2.0)
center = np.array([patch_radius+substrate_thickness, 0, 0])*unit

print("Calculating NF2FF (xz-plane)...")
nf2ff_xz = nf2ff.CalcNF2FF(Sim_Path, f_res, theta, 0, center=center,
                            read_cached=True, outfile='nf2ff_xz.h5')

phi = theta
print("Calculating NF2FF (xy-plane)...")
nf2ff_xy = nf2ff.CalcNF2FF(Sim_Path, f_res, 90, phi, center=center,
                            read_cached=True, outfile='nf2ff_xy.h5')

fig = plt.figure(num="Bent Pattern (polar)", figsize=(15, 7))

ax1 = fig.add_subplot(121, projection='polar')
E_norm_xz = 20.0*np.log10(nf2ff_xz.E_norm/np.max(nf2ff_xz.E_norm)) + nf2ff_xz.Dmax
ax1.plot(np.deg2rad(theta), 10**(np.squeeze(E_norm_xz)/20), linewidth=2, label='xz-plane')
ax1.grid(True)
ax1.set_xlabel('theta (deg)')
ax1.set_theta_zero_location('N')
ax1.set_theta_direction(-1)
ax1.legend(loc=3)

ax2 = fig.add_subplot(122, projection='polar')
E_norm_xy = 20.0*np.log10(nf2ff_xy.E_norm/np.max(nf2ff_xy.E_norm)) + nf2ff_xy.Dmax
ax2.plot(np.deg2rad(phi), 10**(np.squeeze(E_norm_xy)/20), linewidth=2, label='xy-plane')
ax2.grid(True)
ax2.set_xlabel('phi (deg)')
fig.suptitle(f'Bent Patch Antenna Pattern (radius={patch_radius}mm)\nFrequency: {f_res/1e9:.3f} GHz', fontsize=14)
ax2.legend(loc=3)

print(f'Radiated power: Prad = {nf2ff_xy.Prad[0]:.2e} Watt')
print(f'Directivity:    Dmax = {nf2ff_xy.Dmax[0]:.1f} ({10*np.log10(nf2ff_xy.Dmax[0]):.1f} dBi)')
print(f'Efficiency:   nu_rad = {100*nf2ff_xy.Prad[0]/np.real(P_in[idx]):.1f} %')

# --- Full 3D radiation pattern (dome/lobe shape, same style as flat version) ---
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (enables 3D projection)

theta_3d = np.arange(0.0, 180.0, 4.0)
phi_3d   = np.arange(-180.0, 180.0, 4.0)
print("Calculating full 3D far-field (this is the slowest step, be patient)...")
nf2ff_3d = nf2ff.CalcNF2FF(Sim_Path, f_res, theta_3d, phi_3d, center=center)

E_far = nf2ff_3d.E_norm[0] / np.max(nf2ff_3d.E_norm[0])
THETA, PHI = np.meshgrid(theta_3d*np.pi/180, phi_3d*np.pi/180, indexing='ij')

X = E_far * np.sin(THETA) * np.cos(PHI)
Y = E_far * np.sin(THETA) * np.sin(PHI)
Z = E_far * np.cos(THETA)

fig3d = plt.figure(num="3D Radiation Pattern (bent)", figsize=(8, 7))
ax3d = fig3d.add_subplot(111, projection='3d')
ax3d.plot_surface(X, Y, Z, facecolors=plt.cm.jet(E_far),
                   rstride=1, cstride=1, linewidth=0, antialiased=True)
ax3d.set_title(f'3D Radiation Pattern (bent, r={patch_radius}mm) @ {f_res/1e9:.3f} GHz')
ax3d.set_xlabel('X')
ax3d.set_ylabel('Y')
ax3d.set_zlabel('Z')

plt.show()