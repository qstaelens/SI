# Tight-binding hopping integrals for VASP 6.4.2 (and older):
# -----------------------------------------------------------

#TODO: extend for spin-polarized case etc.
#TODO: extend for new VASP 6.5.2 (once available)

import numpy as np
import spglib
from toolbox import *

# nearest neighbours in terms of unit cell (POSCAR): [a1, a2, a3]:
nn_vect = np.array([[0, 0, 0], [0, 1, 0],[0, 0, 1], [1, 0, 0]])
nn_name = np.array([["0"], ["b"], ["c"], ["a"]])

folder = '../'
path_to_POSCAR = folder + 'MLWFS/POSCAR'
path_to_KPOINTS = folder + 'MLWFS/KPOINTS'
path_to_OUTCAR = folder + 'LWL/OUTCAR'
path_to_WANPROJ = folder + 'MLWFS/WANPROJ'
path_to_eigenvalues = folder + 'MLWFS/wannier90.eig'
path_to_INCAR = folder + 'MLWFS/INCAR'

# --- Read system info ---
system, N_b, E_F = read_outcar(path_to_OUTCAR)
print(f'System: {system}')

cell = read_poscar(path_to_POSCAR)
print(f"Space group: {spglib.get_spacegroup(cell, symprec=1e-5)}")

kmesh = read_kpoints(path_to_KPOINTS)
N_k = kmesh[0]*kmesh[1]*kmesh[2]
print(f"Number of kpoints: {N_k}")

excluded, kept_bands = read_exclude_bands(path_to_INCAR, N_b)
print(f"Number of bands: {N_b} --> {len(kept_bands)}")

T, k_points, N_k, N_b, N_w = read_wanproj(path_to_WANPROJ)
print(f"Number of MLWFs: {N_w} \n")

# --- Included bands, ε_nk ---
energies_kept = np.loadtxt(path_to_eigenvalues, usecols=2).reshape(N_k, len(kept_bands))
energies = np.zeros((N_k, N_b))
energies[:, np.array(kept_bands) - 1] = energies_kept

# --- Tightbinding hopping integrals ---
t = np.zeros((len(nn_vect), N_w, N_w))
for R in range(len(nn_vect)):
    for i in range(N_w):
        for j in range(N_w):
            for n in range(N_b):
                for k in range(N_k):
                    t[R, i, j] -= ((np.conj(T[k, n, i])*(energies[k, n]-E_F)*T[k, n, j])*np.exp(1j*2*np.pi*np.dot(k_points[k], nn_vect[R]))).real

    print('Hopping matrix to', nn_name[R], ':')
    with np.printoptions(precision=6, suppress=True, linewidth=1000, floatmode='fixed'):
        matrix_str = np.array2string((1/N_k)*t[R], separator=' ')
        print(matrix_str + '\n')
