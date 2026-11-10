import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colormaps
from matplotlib import rcParams

rcParams.update({'font.size': 20})
rcParams['font.family'] = 'Nimbus Roman'
rcParams["mathtext.fontset"] = "cm"

# Function to extract Fermi energy from OUTCAR
def extract_fermi_energy(file_path):
    with open(file_path, "r") as f:
        for line in f:
            if "E-fermi" in line:
                return float(line.split()[2])  # Third column is the value
    raise ValueError(f"Could not find Fermi energy in {file_path}")

# Function to extract window values from INCAR
def extract_window_values(file_path):
    win_min, win_max = None, None
    with open(file_path, "r") as f:
        for line in f:
            if "dis_win_min" in line:
                win_min = float(line.split("=")[1].strip())
            elif "dis_win_max" in line:
                win_max = float(line.split("=")[1].strip())
    if win_min is None or win_max is None:
        raise ValueError("Could not find window values in INCAR")
    return win_min, win_max

# Read values dynamically
E_F = extract_fermi_energy("../BANDS/OUTCAR")
E_F_MLWF = extract_fermi_energy("../LWL/OUTCAR")
win_min, win_max = extract_window_values("INCAR")

N_bands = 400
N_kseg = 50
K_points = [r'$\Gamma$', 'X', 'M','Z', 'Y', r'$\Gamma$']

bands = True
project = False
windows = False
MLWF = True


if bands: WB = open('../BANDS/OUTCAR', "r")
if project: PC = open('../BANDS/PROCAR', "r")

N_orbs = 9
N_w = 8

# --------------------------------------------------------

WB = WB.read()
WB = np.array([x for x in WB.split(" ") if x != '' and x != '\n'])
k_start = np.where(WB == 'occupation')[0] + 2

E = np.zeros((len(k_start), N_bands))
sum_dk = 0
dk = [0]

for i in range(0,len(k_start)):
    dk += [np.sqrt((float(WB[k_start[i]-9])-float(WB[k_start[i-1]-9]))**2+(float(WB[k_start[i]-8])-float(WB[k_start[i-1]-8]))**2+(float(WB[k_start[i]-7])-float(WB[k_start[i-1]-7]))**2)]
    sum_dk += dk[i]
    e = np.zeros(N_bands)
    for j in range(N_bands):
        e[j] = WB[k_start[i]+3*j]
    E[i] = e
    
dk = (1/sum_dk)*np.array(dk)
E = np.transpose(E)

kk = np.zeros(len(k_start))
for i in range(1,len(k_start)):
    if i == 0:
        kk[0] = 0
    else:
        kk[i] = kk[i-1] + dk[i]
    
K_values = [kk[0]]
for K in range(1, len(K_points)):
    K_values += [kk[K*N_kseg-1]]

if project:
    P = PC.read()
    P = np.array([x for x in P.split(" ") if x != '' and x != '\n'])

    N_k = int(P[5])
    N_b = int(P[9])
    N_i = int(P[13])

    Data = np.zeros((N_b, N_k, N_i, 5))
    for k in range(N_k):
        for b in range(N_b):
            for i in range(N_i):
                # Data[b, k, i, 0] = float(P[32 + 2*N_orbs  + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k].split("\n")[0])
                # Data[b, k, i, 1] = float(P[32 + N_orbs + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])
                for o in [0, 1, 2]:
                    if N_orbs >= 4: Data[b, k, i, 2] += float(P[33 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])
                for o in [0]:
                    if N_orbs >= 9: Data[b, k, i, 3] += float(P[36 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])
                for o in [2]:
                    if N_orbs >= 9: Data[b, k, i, 4] += float(P[36 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])

if MLWF:
    WBI = open('wannier90_band.dat', "r")
    INF = open('wannier90_band.labelinfo.dat', "r")

    WBI = WBI.read()
    WBI = np.array([x for x in WBI.split(" ")[2:] if x != '' and x != '\n'])

    INF = INF.read()
    INF = " ".join(INF.split())
    INF= np.array([x for x in INF.split(" ")])

    k_indices = np.array([int(INF[6*i+1]) for i in range(len(INF)//6)])
    k_coords = np.array([float(INF[6*i+2]) for i in range(len(INF)//6)])

    segment_lengths = np.diff(k_coords)
    total_length = np.sum(segment_lengths)
    num_segments = len(segment_lengths)
    rescaled_coords = [k_coords[0]] 

    for i in range(num_segments):
        frac_tot = kk[N_kseg*(i+1)-1]-kk[N_kseg*(i)]
        segment_points = np.linspace(0, total_length*frac_tot, k_indices[i+1] - k_indices[i])
        rescaled_coords.extend(rescaled_coords[-1] + segment_points)

    rescaled_coords = np.array(rescaled_coords)

    k_full = np.array([float(WBI[2*i]) for i in range(len(WBI)//2)])
    E_full = np.array([float(WBI[2*i+1][:-1]) for i in range(len(WBI)//2)])

    Bands = np.array_split(E_full, len(k_full))
    Bands_full = [float(i[0]) for i in Bands]

    Bands = np.array([Bands_full[i*len(Bands_full)//N_w:i*len(Bands_full)//N_w+len(Bands_full)//N_w] for i in range(N_w)])
    k = rescaled_coords*(K_values[-1]/rescaled_coords[-1])


# Plotting (lm-decomposed) band structures and partial DOS
# --------------------------------------------------------

E -= E_F

fig, axs = plt.subplots(figsize=(8, 6))
norm = plt.Normalize(0, 1)
original_cmap = colormaps.get_cmap('jet')
truncated_cmap = LinearSegmentedColormap.from_list(
    'truncated_jet', original_cmap(np.linspace(0.1, 0.9, 256))
)

if bands:
    for i in range(len(E)):
        # Segment creation
        seg = np.array([[[kk[j], E[i][j]], [kk[j + 1], E[i][j + 1]]] for j in range(len(kk) - 1)])
        
        if project:
            upper_contrib = Data[i, :, 3:5, 3].sum(axis=1)
            lower_contrib = Data[i, :, 3:5, 4].sum(axis=1)
            relative_weight = upper_contrib / (upper_contrib + lower_contrib + 1e-8)            
            total_contrib = upper_contrib + lower_contrib
            whiteness = 1 - total_contrib

        lc = LineCollection(seg, cmap=truncated_cmap)  
        
        if project:
            lc.set_array(relative_weight) 
            lc.set_alpha(1 - whiteness) 
        
        lc.set_linewidth(3.0)
        axs.add_collection(lc)
        line = lc

        plain_lc = LineCollection(seg, colors='black')
        plain_lc.set_linewidth(0.5)
        axs.add_collection(plain_lc)


if MLWF:
    for i in range(len(Bands)):
        if i == 0:
            plt.plot(k[:], Bands[i] - E_F_MLWF, color='black', linewidth=2.0, linestyle='--', label='MLWF')
        else:
            plt.plot(k, Bands[i] - E_F_MLWF, color='black', linewidth=2.0, linestyle='--')
            
if bands and project:
    cbar = fig.colorbar(line, ax=axs, pad=0.08)
    cbar.set_ticks([line.norm.vmin, line.norm.vmax]) 
    cbar.set_ticklabels(['', ''])
    
    # Manually move ticks to top and bottom of the colorbar
    # cbar.ax.tick_params(axis='y', which='both', left=False, right=False) 
    
    cbar.ax.annotate(r'$d_{z^2}$', xy=(0.6, -0.08), xycoords='axes fraction', ha='center', va='bottom')
    cbar.ax.annotate(r'$d_{x^2-y^2}$', xy=(0.6, 1.08), xycoords='axes fraction', ha='center', va='top')


# Plot formatting
plt.ylabel(r'$E-E_{F}$ [eV]')
plt.xlabel('K-points')
plt.ylim(-2, 2)
plt.axhline(0, linestyle=(0, (5, 5)), linewidth=1, color='g', alpha=1)

if windows:
    plt.axhline(win_max - E_F_MLWF, linestyle=(0, (5, 5)), linewidth=1, color='r', alpha=1)
    plt.axhline(win_min - E_F_MLWF, linestyle=(0, (5, 5)), linewidth=1, color='r', alpha=1)

plt.yticks(fontsize=18)
plt.xticks(K_values, K_points)
plt.grid(axis='x')

# Adjust the x-axis to make the plot touch the edges
axs.set_xlim(kk[0], kk[-1])

# Adjust margins and layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.tight_layout(pad=0)

# Save the figure
plt.savefig("bands.pdf", bbox_inches='tight')

