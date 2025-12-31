import numpy as np
import Bio.PDB as PDB
from structure.atomic_protein import ProteinAtomic
from structure.atomic_pigment import ChlorophyllAtomic
from structure.atomic_couplings import calculate_tresp_coupling
from pymembrane.util.physical_constants import kB, hbar, c
from matplotlib import pyplot as plt
from joblib import Parallel, delayed

# Prepare CP43 atomic
# -------------------
path_PSIIcc_pdb = '/Users/doranbennett/Dropbox/Writing2020_2021/Grants/DOE_BSE_PS/2021/Preliminary Data/RengerLineshapes/CrystalStructures/3ARC.pdb'
pdb_PSIIcc = PDB.PDBParser().get_structure('3ARC', path_PSIIcc_pdb)
PSIIcc_atomic = ProteinAtomic(pdb_PSIIcc, name='PSIIcc', )
PSIIcc_atomic.prepare_pigments('CLA', ChlorophyllAtomic, coupling_data_dict='CLA_Renger_2006_jphyschemB_B3LYP_cp43')
cp43_atomic = PSIIcc_atomic.construct_subcomplex(['C'], 'CP43')
cp43_atomic.prepare_pigments('CLA', ChlorophyllAtomic, coupling_data_dict='CLA_Renger_2006_jphyschemB_B3LYP_cp43')

# Prepare Hamiltonian
# -------------------
n_site = len(cp43_atomic.dict_pigments)
site_energy = [14970, 15150, 14790, 14690, 15000, 14940, 14770, 14900, 14770, 14920, 14770, 14860, 14950]
site_label = [f'C_CLA_{index}' for index in np.arange(628, 641)]
site_domain = [0, 1, 1, 1, 2, 3, 2, 2, 2, 4, 2, 5, 5]
list_w_vib = [100, 175, 250, 300, 375, 500, 600, 725, 800, 875]
list_s_vib = [0.2, 0.1, 0.06, 0.04, 0.06, 0.04, 0.015, 0.04, 0.02, 0.02]
H = np.zeros([n_site, n_site], dtype=np.float)
for (index_a, label_a) in enumerate(site_label):
    for (index_b, label_b) in enumerate(site_label):
        if index_a is index_b:
            H[index_a, index_a] = site_energy[index_a]
        else:
            H[index_a, index_b] = calculate_tresp_coupling(cp43_atomic.dict_pigments[label_a],
                                                           cp43_atomic.dict_pigments[label_b],
                                                           0.64, 1.0)

# Double check our domain assignments
dict_domain_labels = {}
for (label, domain) in zip(site_label, site_domain):
    if domain in dict_domain_labels.keys():
        dict_domain_labels[domain].append(label)
    else:
        dict_domain_labels[domain] = [label]
print(dict_domain_labels)


def calculate_00_absorption(index, cp43_atomic, H2_cp43, central_freq):
    # Double check our domain assignments
    dict_domain_labels = {}
    for (label, domain) in zip(site_label, site_domain):
        if domain in dict_domain_labels.keys():
            dict_domain_labels[domain].append(label)
        else:
            dict_domain_labels[domain] = [label]
    print(dict_domain_labels)

    # Build list of sites for each domain
    dict_domain_index = {}
    for domain in set(site_domain):
        dict_domain_index[domain] = [site_label.index(label) for label in dict_domain_labels[domain]]

    # Construct the domain coefficients
    # ---------------------------------
    dict_domain_coefficients = {}
    dict_domain_energies = {}
    for domain in set(site_domain):
        list_sites = dict_domain_index[domain]
        output = np.linalg.eigh(H2_cp43[np.ix_(list_sites, list_sites)])
        # coefficients are listed such that the column is the eigenvector of H
        dict_domain_coefficients[domain] = output[1]
        dict_domain_energies[domain] = output[0]

    # Prepare Renger Site Lineshape Functions
    # -----------------------------------------
    dt = 0.1
    t_max = 2000.0
    t_axis = np.arange(0, t_max, dt)

    dw = 0.1
    w_max = 2000
    w_axis = np.arange(dw, w_max, dw)

    # Spectral density
    def J_w(w_axis):
        S0 = 0.5
        s1 = 0.8
        s2 = 0.5
        w1 = 0.069 * 8.0656  # Units: cm^-1
        w2 = 0.24 * 8.0656  # Units: cm^-1
        prefactor = S0 / (s1 + s2) * np.power(w_axis, 3) / 10080
        sd_1 = s1 / w1 ** 4 * np.exp(-np.sqrt(w_axis / w1))
        sd_2 = s2 / w2 ** 4 * np.exp(-np.sqrt(w_axis / w2))
        return prefactor * (sd_1 + sd_2)

    Jw_renger = J_w(w_axis)
    temp = 77.0  # Units: K
    e_lambda = np.trapz(w_axis * Jw_renger, x=w_axis)

    def n_w(w_axis, temp):
        return 1 / (np.exp(w_axis / (kB * temp)) - 1)

    def g_t(t_axis, w_axis, j_w, temp):
        return np.array([np.trapz((1 + n_w(w_axis, temp)) * j_w * np.exp(-1j * w_axis / hbar * t)
                                  + n_w(w_axis, temp) * j_w * np.exp(1j * w_axis / hbar * t),
                                  x=w_axis) for t in t_axis])

    g_t = g_t(t_axis, w_axis, Jw_renger, temp)

    # Calculate Lifetime broadening
    # -----------------------------
    # I have set the correlation distance to 0 in order to simplify the calculation
    # of lifetimes.
    dict_domain_lifetimes = {}
    for domain in set(site_domain):
        list_tau = []
        for exciton in np.arange(len(dict_domain_energies[domain])):
            tau_exciton = 0
            for exciton_N in np.arange(len(dict_domain_energies[domain])):
                if exciton_N != exciton:
                    omega_MN = dict_domain_energies[domain][exciton] - dict_domain_energies[domain][exciton_N]
                    gamma_MN = np.dot(np.abs(dict_domain_coefficients[domain][:, exciton]) ** 2,
                                      np.abs(dict_domain_coefficients[domain][:, exciton_N]) ** 2)
                    if omega_MN > 0:
                        tau_exciton += np.pi * gamma_MN * omega_MN ** 2 * (1 + n_w(omega_MN, temp)) * J_w(omega_MN)
                    else:
                        tau_exciton += np.pi * gamma_MN * omega_MN ** 2 * (n_w(-omega_MN, temp)) * J_w(-omega_MN)
            list_tau.append(tau_exciton)
        dict_domain_lifetimes[domain] = list_tau

    # Construct exciton lineshape components
    # --------------------------------------
    dict_exciton_gamma = {}
    for domain in set(site_domain):
        dict_exciton_gamma[domain] = np.sum(np.abs(dict_domain_coefficients[domain]) ** 4, axis=0)

    dict_exciton_etrans = {}
    for domain in set(site_domain):
        dict_exciton_etrans[domain] = [
            dict_domain_energies[domain][exciton] - dict_exciton_gamma[domain][exciton] * e_lambda
            for exciton in np.arange(len(dict_exciton_gamma[domain]))]

    # Construct Absorption Spectrum for each exciton
    # ----------------------------------------------
    dict_exciton_Dw = {}
    for domain in set(site_domain):
        list_Dw_by_exciton = []
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            w_Md = dict_exciton_etrans[domain][exciton] - central_freq
            G_Md = dict_exciton_gamma[domain][exciton] * g_t
            Dt = np.exp(-1j * w_Md * t_axis / hbar) * np.exp(G_Md - G_Md[0]) * np.exp(
                -t_axis * dict_domain_lifetimes[domain][exciton] / hbar)
            list_Dw_by_exciton.append(1 / (2 * np.pi) * np.fft.fft(Dt))

        dict_exciton_Dw[domain] = list_Dw_by_exciton

    # Calculate the exciton dipole moment
    # -----------------------------------
    dict_domain_dipoles = {}
    for domain in set(site_domain):
        list_dipole_by_exciton = []
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            mu_exc = np.zeros(3)
            for (index_site, label_site) in enumerate(dict_domain_labels[domain]):
                mu_exc += dict_domain_coefficients[domain][index_site, exciton] * cp43_atomic.dict_pigments[
                    label_site].get_dipole_dir()

            list_dipole_by_exciton.append(mu_exc)

        dict_domain_dipoles[domain] = list_dipole_by_exciton

    # Calculate Absorption Spectrum
    # -----------------------------
    A_w = np.zeros(len(t_axis))
    for domain in set(site_domain):
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            A_w += np.linalg.norm(dict_domain_dipoles[domain][exciton]) ** 2 * np.real(dict_exciton_Dw[domain][exciton])

    path_save = '/Users/doranbennett/Dropbox/Writing2020_2021/Grants/DOE_BSE_PS/2021/Preliminary Data/RengerLineshapes/LineshapeTests/CP43_2015'
    np.save(path_save + f'/Aw_00_{index}.npy', A_w)


def calculate_01_absorption(index, cp43_atomic, H2_cp43, central_freq):
    # Prepare Renger Site Lineshape Functions
    # -----------------------------------------
    dt = 0.1
    t_max = 2000.0
    t_axis = np.arange(0, t_max, dt)

    dw = 0.1
    w_max = 2000
    w_axis = np.arange(dw, w_max, dw)

    # Spectral density
    def J_w(w_axis):
        S0 = 0.5
        s1 = 0.8
        s2 = 0.5
        w1 = 0.069 * 8.0656  # Units: cm^-1
        w2 = 0.24 * 8.0656  # Units: cm^-1
        prefactor = S0 / (s1 + s2) * np.power(w_axis, 3) / 10080
        sd_1 = s1 / w1 ** 4 * np.exp(-np.sqrt(w_axis / w1))
        sd_2 = s2 / w2 ** 4 * np.exp(-np.sqrt(w_axis / w2))
        return prefactor * (sd_1 + sd_2)

    Jw_renger = J_w(w_axis)
    temp = 77.0  # Units: K
    e_lambda = np.trapz(w_axis * Jw_renger, x=w_axis)

    def n_w(w_axis, temp):
        return 1 / (np.exp(w_axis / (kB * temp)) - 1)

    def g_t(t_axis, w_axis, j_w, temp):
        return np.array([np.trapz((1 + n_w(w_axis, temp)) * j_w * np.exp(-1j * w_axis / hbar * t)
                                  + n_w(w_axis, temp) * j_w * np.exp(1j * w_axis / hbar * t),
                                  x=w_axis) for t in t_axis])

    g_t = g_t(t_axis, w_axis, Jw_renger, temp)

    list_Dw_site = []
    Aw_01 = 0j*t_axis
    site_energy = np.diag(H2_cp43)
    for site_index in np.arange(len(cp43_atomic.dict_pigments)):
        Dw_site = 0*g_t
        for (w_vib, s_vib) in zip(list_w_vib, list_s_vib):
            fc_01_sq = np.exp(-s_vib)*(s_vib)
            fc_00_sq = np.exp(-s_vib)
            w_mi = site_energy[site_index]+w_vib-e_lambda - central_freq
            Dw_site += fc_01_sq/(fc_00_sq*2*np.pi)*np.real(np.fft.fft(np.exp(-1j*w_mi*t_axis/hbar)*np.exp(g_t - g_t[0])))
        list_Dw_site.append(Dw_site)
        Aw_01 += Dw_site

    path_save = '/Users/doranbennett/Dropbox/Writing2020_2021/Grants/DOE_BSE_PS/2021/Preliminary Data/RengerLineshapes/LineshapeTests/CP43_2015'
    np.save(path_save + f'/Aw_01_{index}.npy', Aw_01)
    print(f'Index {index} finished vibrational contributions. ')


# Construct the Absorption Spectrum
# =================================
N_ens = 100
sigma_e = 120/(2*np.sqrt(2*np.log(2)))

def calculate_spectra_w_vib(index, cp43_atomic, H2_cp43):
    rng = np.random.RandomState(seed=index)
    H_wdis = H2_cp43 + np.diag(rng.normal(0,
                                 scale=sigma_e,
                                 size=len(site_energy)))
    central_freq = np.mean(site_energy)
    calculate_00_absorption(index, cp43_atomic, H_wdis, central_freq)
    calculate_01_absorption(index, cp43_atomic, H_wdis, central_freq)


Parallel(n_jobs=-3)(delayed(calculate_spectra_w_vib)(n, cp43_atomic,H) for n in range(N_ens))

for index in range(N_ens): 
    Aw_00 = np.load(
        f'/Users/doranbennett/Dropbox/Writing2020_2021/Grants/DOE_BSE_PS/2021/Preliminary Data/RengerLineshapes/LineshapeTests/CP43_2015/Aw_00_{index}.npy')
    Aw_01 = np.load(
        f'/Users/doranbennett/Dropbox/Writing2020_2021/Grants/DOE_BSE_PS/2021/Preliminary Data/RengerLineshapes/LineshapeTests/CP43_2015/Aw_01_{index}.npy')
    if index == 0:
        Aw_mean = np.real(Aw_00)+np.real(Aw_01)
    else:
        Aw_mean += np.real(Aw_00)+np.real(Aw_01)

# Construct frequency axis
central_freq = np.mean(site_energy)
dt = 0.1
t_max = 2000.0
t_axis = np.arange(0, t_max, dt)
w_axis = -(2*np.pi*hbar*np.fft.fftfreq(len(t_axis), dt)-central_freq)
lambda_axis = c/(w_axis/(2*np.pi*hbar))


plt.plot(np.fft.fftshift(lambda_axis), np.fft.fftshift((Aw_mean) / np.max(Aw_mean)), 'k')
plt.xlim(650, 700)