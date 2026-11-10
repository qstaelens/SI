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


# Read values dynamically
E_F = extract_fermi_energy("../BANDS/OUTCAR")

N_bands = 400
N_kseg = 50
K_points = [r'$\Gamma$', 'X', 'M','Z', 'Y', r'$\Gamma$']

bands = True
project = True


if bands: WB = open('../BANDS/OUTCAR', "r")
if project: PC = open('../BANDS/PROCAR', "r")

N_orbs = 9

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

    Data = np.zeros((N_b, N_k, N_i, 6))
    for k in range(N_k):
        for b in range(N_b):
            for i in range(N_i):
                Data[b, k, i, 0] = float(P[32 + 2*N_orbs  + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k].split("\n")[0])   #total
                Data[b, k, i, 1] = float(P[32 + N_orbs + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k])                     #s orbital
                for o in [0, 1, 2]:
                    if N_orbs >= 4: Data[b, k, i, 2] += float(P[33 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k]) #all p orbitals
                for o in [2]:
                    if N_orbs >= 4: Data[b, k, i, 3] += float(P[33 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k]) #px orbital
                for o in [0, 1, 2, 3, 4]:
                    if N_orbs >= 9: Data[b, k, i, 4] += float(P[36 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k]) #all d orbitals
                for o in [4]:
                    if N_orbs >= 9: Data[b, k, i, 5] += float(P[36 + N_orbs + o + (N_orbs + 2)*i + ((N_i + 2)*(N_orbs + 2) + 5)*b + (((N_i + 2)*(N_orbs + 2) + 5)*N_b+ 9)*k]) #dx2-y2 orbital


print(Data)

# Plotting (lm-decomposed) band structures and partial DOS
# --------------------------------------------------------

E -= E_F

fig, axs = plt.subplots(figsize=(8, 6))
norm = plt.Normalize(0, 1)
cmap = plt.get_cmap("Blues")
full_cmap = plt.get_cmap("Blues")
truncated_cmap = LinearSegmentedColormap.from_list(
    'truncated_Blues', full_cmap(np.linspace(0.2, 1, 256))
)

full_cmap = plt.get_cmap("Reds")
truncated_cmap = LinearSegmentedColormap.from_list(
    'truncated_Reds', full_cmap(np.linspace(0.2, 1, 256))
)

if bands:
    for i in range(len(E)):
        # Segment creation
        seg = np.array([[[kk[j], E[i][j]], [kk[j + 1], E[i][j + 1]]] for j in range(len(kk) - 1)])
        
        if project:
            p = Data[i, :, 24:52, 2].sum(axis=1)
            px = Data[i, :, 24:52, 3].sum(axis=1)
            d = Data[i, :, 16:24, 4].sum(axis=1)
            dx2y2 = Data[i, :, 16:24, 5].sum(axis=1)
            
            q = dx2y2

        lc = LineCollection(seg, cmap=full_cmap, norm=norm)  
        
        if project:
            lc.set_array(q) 
            lc.set_alpha(q) 
        
        lc.set_linewidth(3.0)
        axs.add_collection(lc)
        line = lc

        plain_lc = LineCollection(seg, colors='black')
        plain_lc.set_linewidth(0.3)
        axs.add_collection(plain_lc)


if bands and project:
    cbar = fig.colorbar(line, ax=axs, pad=0.08)
    ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])


# Plot formatting
plt.ylabel(r'$E-E_{F}$ [eV]')
plt.xlabel('K-points')
plt.ylim(-2.0, 2.0)
plt.axhline(0, linestyle=(0, (5, 5)), linewidth=1, color='g', alpha=1)
plt.title("Cu: $d_{x^2-y^2}$", fontsize=20, pad=20)
#plt.title("Cu: $d$", fontsize=20, pad=20)
#plt.title("O: $p_{x}$", fontsize=20, pad=20)
#plt.title("O: $p$", fontsize=20, pad=20)

plt.yticks(fontsize=18)
plt.xticks(K_values, K_points)
plt.grid()

# Adjust the x-axis to make the plot touch the edges
axs.set_xlim(kk[0], kk[-1])

# Adjust margins and layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.tight_layout(pad=0)
plt.savefig("LM-dx2y2.png", bbox_inches='tight')

