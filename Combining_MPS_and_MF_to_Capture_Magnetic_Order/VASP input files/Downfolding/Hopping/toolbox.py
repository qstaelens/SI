#TODO: clean all and remove redundancies

import numpy as np
import random

def read_poscar(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    scaling_factor = float(lines[1].strip())

    # Lattice vectors
    lattice = []
    for i in range(2, 5):
        lattice.append([float(x) for x in lines[i].strip().split()])
    lattice = np.array(lattice) * scaling_factor

    # Element counts
    element_counts = [int(x) for x in lines[6].strip().split()]
    total_atoms = sum(element_counts)

    # Atomic numbers by type (1, 2, 3, ...)
    numbers = []
    atom_type = 1
    for count in element_counts:
        numbers.extend([atom_type] * count)
        atom_type += 1

    # Positions (start reading from line 8, which is line index 7)
    fractional_positions = []
    position_start = 8
    for line in lines[position_start:position_start + total_atoms]:
        tokens = line.strip().split()
        if len(tokens) >= 3:
            fractional_positions.append([float(x) for x in tokens[:3]])

    return lattice, fractional_positions, numbers

def read_kpoints(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    mesh = [int(x) for x in lines[3].strip().split()]
    
    return mesh

def read_outcar(filename):
    system = None
    N_b = None
    E_F = None

    with open(filename, 'r') as file:
        for line in file:
            if 'SYSTEM =' in line:
                system = line.split('=')[1].strip()
            elif 'NBANDS =' in line:
                N_b = int(line.split('=')[1].strip())
            elif 'E-fermi' in line:  # Fermi level is usually written like: "E-fermi :   5.1234"
                try:
                    E_F = float(line.split()[2])
                except (IndexError, ValueError):
                    pass  # skip if parsing fails

    return system, N_b, E_F


def read_wanproj(path_to_WANPROJ: str):
    """
    Reads a VASP WANPROJ file and extracts the Wannier transformation matrices (T)
    and k-points, automatically reading N_k, N_b, N_w from the file.
    
    Parameters:
        path_to_WANPROJ (str): Path to the WANPROJ file.
    
    Returns:
        T (np.ndarray): Wannier transformation matrix of shape (N_k, N_b, N_w), complex dtype.
        k_points (np.ndarray): Array of k-points of shape (N_k, 3), float dtype.
        N_k (int): Number of k-points.
        N_b (int): Number of Bloch bands.
        N_w (int): Number of Wannier functions.
    """
    
    # Read all lines
    with open(path_to_WANPROJ, "r") as f:
        lines = f.readlines()
    
    # Skip comment lines
    data_lines = [line for line in lines if not line.strip().startswith("#")]
    
    # First non-comment line contains 1 N_k N_b N_w
    first_numbers = data_lines[0].split()
    N_k = int(first_numbers[1])
    N_b = int(first_numbers[2])
    N_w = int(first_numbers[3])
    
    # Flatten all remaining data into a single list of strings
    proj_data = []
    for line in data_lines[1:]:
        proj_data.extend(line.split())  # split handles spaces/tabs
    
    # Extract k-points (next 4*N_k entries)
    k_data = proj_data[:4 * N_k]
    k_points = np.array([float(k_data[i]) for i in range(len(k_data))]).reshape(N_k, 4)[:, 1:]
    
    # Remaining data for T
    T_data = proj_data[4 * N_k:]
    
    T = np.zeros((N_k, N_b, N_w), dtype=complex)
    N_bk_index = 1
    
    for K in range(N_k):
        N_bk = int(T_data[N_bk_index])
        for b in range(N_bk):
            N_b_index = int(T_data[4 + 4 * N_w * b + N_bk_index]) - 1
            for i in range(N_w):
                Re = float(T_data[6 + 4 * N_w * b + 4 * i + N_bk_index])
                Im = float(T_data[7 + 4 * N_w * b + 4 * i + N_bk_index])
                T[K, N_b_index, i] = complex(Re, Im)
        N_bk_index += 4 * N_w * N_bk + 5
    
    return T, k_points, N_k, N_b, N_w



def read_exclude_bands(filename, N_b):

    exclude_bands_line = None

    with open(filename, "r") as f:
        for line in f:
            if line.strip().lower().startswith("exclude_bands"):
                exclude_bands_line = line.split("=")[1].strip()
                break

    excluded = set()
    if exclude_bands_line:
        ranges = exclude_bands_line.split(",")
        for r in ranges:
            if "-" in r:
                start, end = map(int, r.split("-"))
                excluded.update(range(start, end + 1))
            else:
                excluded.add(int(r))

    all_bands = np.arange(1, N_b + 1)
    kept_bands = [b for b in all_bands if b not in excluded]

    return excluded, kept_bands


# Function to convert fractional to Cartesian coordinates
def fractional_to_cartesian(fractional_coord, a_vec, b_vec, c_vec):
    return fractional_coord[0] * a_vec + fractional_coord[1] * b_vec + fractional_coord[2] * c_vec


def Perturb_my_POSCAR(file_path, output_file_path, max_shift):
    modified_lines = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        in_direct_coords = False  # Flag to check if we are in the Direct coordinate section

        for line in lines:
            # Detect the start of the Direct coordinates
            if line.strip().lower().startswith("direct"):
                in_direct_coords = True
                modified_lines.append(line)
                continue

            # Check if we are in the Direct coordinates section and if the line contains coordinates
            if in_direct_coords:
                # Skip empty lines or non-coordinate lines
                if line.strip() == "":
                    modified_lines.append(line)
                    continue

                # Parse the line and add a random float to each coordinate value
                coords = line.split()
                if len(coords) >= 3:  # Expecting at least three values for x, y, and z coordinates
                    modified_coords = [
                        f"{float(coord) + random.uniform(0.0, max_shift):.16f}" for coord in coords[:3]
                    ]
                    modified_line = "  ".join(modified_coords) + "\n"
                    modified_lines.append(modified_line)
                else:
                    # If line does not contain 3 coordinates, append it unchanged
                    modified_lines.append(line)
            else:
                # Before reaching Direct coordinates, add lines unmodified
                modified_lines.append(line)

    # Write the modified lines to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(modified_lines)



import numpy as np

def read_poscar_lattice(file_path):
    """Reads the lattice vectors from a POSCAR file."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Get the scaling factor from line 2
    scale = float(lines[1].strip())
    
    # Read lattice vectors from lines 3 to 5 and apply scaling
    lattice_vectors = []
    for i in range(3):
        vector = np.array([float(x) for x in lines[2 + i].strip().split()])
        lattice_vectors.append(vector * scale)
    
    # Convert lattice_vectors list to a numpy array
    lattice_matrix = np.array(lattice_vectors)
    return lattice_matrix

def eigsolve_lattice_vectors(lattice_matrix):
    """Computes the eigenvalues and eigenvectors of the lattice matrix."""
    eigenvalues, eigenvectors = np.linalg.eig(lattice_matrix)
    return eigenvalues, eigenvectors
