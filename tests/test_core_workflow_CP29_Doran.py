import numpy as np
import pandas as pd
import Bio.PDB as PDB
from pymembrane.structure.atomic_protein import ProteinAtomic
from pymembrane.structure.atomic_pigment import ChlorophyllAtomic
from pymembrane.structure.atomic_couplings import calculate_tresp_coupling
from pymembrane.util.physical_constants import kB, hbar, c
from matplotlib import pyplot as plt

# Inputs
# ------
dt = 0.1
t_max = 2000.0
t_axis = np.arange(0, t_max, dt)

dw = 0.1
w_max = 2000
w_axis = np.arange(dw, w_max, dw)

temp = 130.0  # Units: K

N_ens = 300
sigma_e = 130/(2*np.sqrt(2*np.log(2)))

# Prepare CP29 atomic
# -------------------
working_dir = '/Users/doranbennett/Dropbox/Writing2022/Manuscripts/PSII_Supercomplex/CP29RengerSpectra/'
pdb_name = working_dir + 'pqr_extended.pdb'
cp29_atomic = ProteinAtomic(f'{pdb_name}', "cp29_monomer", center_atomic=False,
                            set_atomic_charges=True)
cp29_atomic.prepare_pigments('CLA', ChlorophyllAtomic,
                             coupling_data_dict='CLA_Kleinekathofer_2021_JCP_B3LYP_cp29')
cp29_atomic.prepare_pigments('CHL', ChlorophyllAtomic,
                             coupling_data_dict='CHL_Kleinekathofer_2021_JCP_B3LYP_cp29')

# Prepare Hamiltonian
# -------------------
n_site = len(cp29_atomic.dict_pigments)

# e0 = 14920

# fitted site energy form Renger's paper:
# =======================================
site_energy = [14980, 14900, 14900, 15684, 15439, 15439, 14980, 14920, 14850, 14900,
               14880, 15674, 14940]
central_freq = np.mean(site_energy)

list_site_label = ['C_CLA_602', 'C_CLA_603', 'C_CLA_604', 'C_CHL_606',
                   'C_CHL_607', 'C_CHL_608', 'C_CLA_609', 'C_CLA_610',
                   'C_CLA_611', 'C_CLA_612', 'C_CLA_613', 'C_CHL_614',
                   'C_CLA_615']
list_site_domain = [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 2, 3, 0]
list_w_vib = [100, 175, 250, 300, 375, 500, 600, 725, 800, 875]
list_s_vib = [0.2, 0.1, 0.06, 0.04, 0.06, 0.04, 0.015, 0.04, 0.02, 0.02]

# H = np.zeros([n_site, n_site], dtype=np.float64)
H = pd.read_csv(f"{working_dir}H_cp29_Renger_fitted_sEnergy.csv")
H = np.array(H.values)

# Double check our domain assignments
mu_cla = 5.47 # Units: Debye
mu_chl = 4.61 # Units: Debye
dict_dipole_by_site = {'C_CLA_602':mu_cla,
                      'C_CLA_603':mu_cla,
                      'C_CLA_604':mu_cla,
                      'C_CHL_606':mu_chl,
                      'C_CHL_607':mu_chl,
                      'C_CHL_608':mu_chl,
                      'C_CLA_609':mu_cla,
                      'C_CLA_610':mu_cla,
                      'C_CLA_611':mu_cla,
                      'C_CLA_612':mu_cla,
                      'C_CLA_613':mu_cla,
                      'C_CHL_614':mu_chl,
                      'C_CLA_615':mu_cla}

# Double check our domain assignments
dict_domain_labels = {}
for (label, domain) in zip(list_site_label, list_site_domain):
    if domain in dict_domain_labels.keys():
        dict_domain_labels[domain].append(label)
    else:
        dict_domain_labels[domain] = [label]
print(dict_domain_labels)

# Prepare Renger Site Lineshape Functions
# -----------------------------------------

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
e_lambda = np.trapz(w_axis * Jw_renger, x=w_axis)
print('e_lambda is ', e_lambda)


def n_w(w_axis, temp):
    return 1 / (np.exp(w_axis / (kB * temp)) - 1)


def g_t(t_axis, w_axis, j_w, temp):
    return np.array(
        [np.trapz((1 + n_w(w_axis, temp)) * j_w * np.exp(-1j * w_axis / hbar * t)
                  + n_w(w_axis, temp) * j_w * np.exp(1j * w_axis / hbar * t),
                  x=w_axis) for t in t_axis])


g_t = g_t(t_axis, w_axis, Jw_renger, temp)

def calculate_spectra(H2, atomic_protein):
    # Build list of sites for each domain
    dict_domain_index = {}
    for domain in set(list_site_domain):
        dict_domain_index[domain] = [list_site_label.index(label) for label in
                                     dict_domain_labels[domain]]

    # Construct the Domain Coefficients
    # ---------------------------------
    dict_domain_coefficients = {}
    dict_domain_energies = {}
    for domain in set(list_site_domain):
        list_sites = dict_domain_index[domain]
        output = np.linalg.eigh(H2[np.ix_(list_sites, list_sites)])
        # coefficients are listed such that the column is the eigenvector of H
        dict_domain_coefficients[domain] = output[1]
        dict_domain_energies[domain] = output[0]

    # Calculate Lifetime Broadening
    # -----------------------------
    # I have set the correlation distance to 0 in order to simplify the calculation
    # of lifetimes.
    dict_domain_lifetimes = {}
    for domain in set(list_site_domain):
        list_tau = []
        for exciton in np.arange(len(dict_domain_energies[domain])):
            tau_exciton = 0
            for exciton_N in np.arange(len(dict_domain_energies[domain])):
                if exciton_N != exciton:
                    omega_MN = dict_domain_energies[domain][exciton] - \
                               dict_domain_energies[domain][exciton_N]
                    gamma_MN = np.dot(
                        np.abs(dict_domain_coefficients[domain][:, exciton]) ** 2,
                        np.abs(dict_domain_coefficients[domain][:, exciton_N]) ** 2)
                    if omega_MN > 0:
                        tau_exciton += np.pi * gamma_MN * omega_MN ** 2 * (
                                1 + n_w(omega_MN, temp)) * J_w(omega_MN)
                    else:
                        tau_exciton += np.pi * gamma_MN * omega_MN ** 2 * (
                            n_w(-omega_MN, temp)) * J_w(-omega_MN)
            list_tau.append(tau_exciton)
        dict_domain_lifetimes[domain] = list_tau

    # Construct Exciton Lineshape Components
    # --------------------------------------
    dict_exciton_gamma = {}
    for domain in set(list_site_domain):
        dict_exciton_gamma[domain] = np.sum(
            np.abs(dict_domain_coefficients[domain]) ** 4, axis=0)

    dict_exciton_etrans = {}
    for domain in set(list_site_domain):
        dict_exciton_etrans[domain] = [
            dict_domain_energies[domain][exciton] - (
                    dict_exciton_gamma[domain][exciton] * e_lambda)
            for exciton in np.arange(len(dict_exciton_gamma[domain]))]

    # Calculate Boltzmann Thermalization
    # ----------------------------------
    bolt_denom = 0.
    for domain in set(list_site_domain):
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            bolt_denom += np.exp((-dict_exciton_etrans[domain][exciton]) / (kB * temp))

    dict_exciton_therm = {}
    for domain in set(list_site_domain):
        dict_exciton_therm[domain] = [np.exp((-dict_exciton_etrans[domain][exciton]) /
                                             (kB * temp)) / bolt_denom for exciton in
                                      np.arange(len(dict_exciton_gamma[domain]))]

    # Calculate the Exciton Dipole Moment
    # -----------------------------------
    # Note to fix this for different pigment types
    dict_domain_dipoles = {}
    for domain in set(list_site_domain):
        list_dipole_by_exciton = []
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            mu_exc = np.zeros(3)
            for (index_site, label_site) in enumerate(dict_domain_labels[domain]):
                mu_exc += (dict_domain_coefficients[domain][index_site, exciton]
                           * atomic_protein.dict_pigments[
                              label_site].get_dipole_dir()
                           * dict_dipole_by_site[label_site])
            list_dipole_by_exciton.append(mu_exc)
        dict_domain_dipoles[domain] = list_dipole_by_exciton

    # Construct Spectrum for each Exciton
    # -----------------------------------
    dict_exciton_fl_contribution = {}
    dict_exciton_a_contribtion = {}
    d_ti_pre = 0 * g_t
    for (w_vib, s_vib) in zip(list_w_vib, list_s_vib):
        fc_01_sq = np.exp(-s_vib) * (s_vib)
        fc_00_sq = np.exp(-s_vib)
        w_mi = - w_vib - e_lambda - central_freq
        g_mi = g_t - g_t[0]
        d_ti_pre += np.exp(1j * w_mi * t_axis / hbar) * np.exp(g_mi) * fc_01_sq / (fc_00_sq)

    Fw_site = 0 * g_t
    for domain in set(list_site_domain):
        list_fl_contribution_by_exciton = []
        list_a_contribution_by_exciton = []
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            w_Md = dict_exciton_etrans[domain][exciton] - central_freq
            G_Md = dict_exciton_gamma[domain][exciton] * g_t
            # Absorption
            A_cont = 0.
            Dt_a = np.exp(-1j * w_Md * t_axis / hbar) * np.exp(G_Md - G_Md[0]) * np.exp(
                -t_axis * dict_domain_lifetimes[domain][exciton] / hbar)
            Dw_a = 1 / (2 * np.pi) * np.fft.fft(Dt_a)
            A_cont += np.linalg.norm(
                dict_domain_dipoles[domain][exciton]) ** 2 * np.real(Dw_a)
            # Fluorescence
            F_cont = 0.
            Dt_f = np.exp(1j * w_Md * t_axis / hbar) * np.exp(G_Md - G_Md[0]) * np.exp(
                -t_axis * dict_domain_lifetimes[domain][exciton] / hbar)
            Dw_f = 1 / (2 * np.pi) * np.fft.fft(Dt_f)
            F_cont += dict_exciton_therm[domain][exciton] * \
                      np.linalg.norm(
                          dict_domain_dipoles[domain][exciton]) ** 2 * np.real(Dw_f)

            site_energy = np.diag(H2)
            Ft_cont = 0 * g_t
            for site_index in dict_domain_index[domain]:
                site_label = list_site_label[site_index]
                coef_index = \
                    np.where(np.array(dict_domain_index[domain]) == site_index)[0][0]
                coef_sq = dict_domain_coefficients[domain][:, exciton][coef_index] ** 2
                d_ti_site = (dict_exciton_therm[domain][exciton]
                             * dict_dipole_by_site[site_label]**2
                             * coef_sq
                             * np.exp(1j * site_energy[site_index] * t_axis / hbar)
                             * d_ti_pre)
                Ft_cont += d_ti_site
            Fw_site += np.real(np.fft.fft(Ft_cont)/(2 * np.pi))
            list_fl_contribution_by_exciton.append(F_cont)
            list_a_contribution_by_exciton.append(A_cont)
        dict_exciton_fl_contribution[domain] = list_fl_contribution_by_exciton
        dict_exciton_a_contribtion[domain] = list_a_contribution_by_exciton

    Dt_pre = 0 * g_t
    for (w_vib, s_vib) in zip(list_w_vib, list_s_vib):
        fc_01_sq = np.exp(-s_vib) * (s_vib)
        fc_00_sq = np.exp(-s_vib)
        w_mi = w_vib - e_lambda - central_freq
        Dt_pre += fc_01_sq / (fc_00_sq) * np.exp(-1j * w_mi * t_axis / hbar) * np.exp(g_t - g_t[0])

    site_energy = np.diag(H2)
    Dt_site = 0 * g_t
    for site_index in np.arange(len(atomic_protein.dict_pigments)):
        site_label = list_site_label[site_index]
        Dt_site += (dict_dipole_by_site[site_label]**2
                   * np.exp(-1j * site_energy[site_index] * t_axis / hbar)
                   * Dt_pre)
    Aw_site = np.real(np.fft.fft(Dt_site)/(2 * np.pi))

    # Calculate Fluorescence Spectrum
    # -----------------------------
    F_w = Fw_site
    for domain in set(list_site_domain):
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            F_w += dict_exciton_fl_contribution[domain][exciton]

    # Calculate Absorption Spectrum
    # -----------------------------
    A_w = Aw_site
    for domain in set(list_site_domain):
        for exciton in np.arange(len(dict_exciton_gamma[domain])):
            A_w += dict_exciton_a_contribtion[domain][exciton]
    return F_w, A_w, Fw_site, dict_exciton_fl_contribution, Aw_site, dict_exciton_a_contribtion

F_w_vib = 0.
A_w_vib = 0.
for index in np.arange(N_ens):
    rng = np.random.RandomState(seed=index)
    H_wdis = H + np.diag(rng.normal(0,
                                 scale=sigma_e,
                                 size=len(site_energy)))
    central_freq = np.mean(site_energy)
    output = calculate_spectra(H_wdis, cp29_atomic)
    F_w_vib += output[0]
    A_w_vib += output[1]
Fw_norm = F_w_vib/np.max(F_w_vib)
Aw_norm = A_w_vib/np.max(A_w_vib)


########################################
####    Plot Fluorescence Spectra    ###
########################################
spectra_data = '/Users/doranbennett/Dropbox/Writing2022/Manuscripts/PSII_Supercomplex/CP29RengerSpectra/abso_fluor_data/'

# Construct frequency axis
w_axis_f = (2 * np.pi * hbar * np.fft.fftfreq(len(t_axis), dt) + central_freq)
lambda_axis_f = c / (w_axis_f / (2 * np.pi * hbar))

# Experimental Spectra
fluorescence_exp_125k = np.load(spectra_data + 'fluorescence_exp_125k.npy')
fluorescence_renger_125 = np.load(spectra_data + 'fluorescence_renger_125.npy')


fig, (ax) = plt.subplots(1, 1)
ax.plot(np.fft.fftshift(lambda_axis_f)[9500:10100],
        np.fft.fftshift(Fw_norm)[9500:10100], label="Renger's Hamiltonian",
        linewidth=4.5, color='Blue', alpha=0.6)
ax.scatter(fluorescence_exp_125k[:, 0], fluorescence_exp_125k[:, 1] / np.max(
    fluorescence_exp_125k[:, 1]), color='slategray', label='Experiment', linewidth=2,
           linestyle='dashed')
ax.scatter(fluorescence_renger_125[:, 0], fluorescence_renger_125[:, 1] / np.max(
    fluorescence_renger_125[:, 1]), color='green', label='Renger', linewidth=2)
ax.set(xlim=(660, 730))
ax.set_xlabel('Wavelength (nm)', fontsize=18)
ax.set_ylabel(f'Fluorescence (a.u.) at {temp}K, \n N_ens= {N_ens},sigma={sigma_e}',
              fontsize=18)
ax.set_yticks([0, 0.5, 1])
ax.set_xticks([660, 680, 700, 720, 750])
ax.set_yticklabels(labels=[0, 0.5, 1], fontsize=15)
ax.set_xticklabels(labels=[660, 680, 700, 720, 750], fontsize=15)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.legend()
plt.tight_layout()
# plt.savefig(f"{working_dir}flor_{temp}K_N_ens_{N_ens}_sigma_{sigma_e}.png")
fig.show()

########################################
####     Plot Absorption Spectra     ###
########################################

w_axis_a = -(2 * np.pi * hbar * np.fft.fftfreq(len(t_axis), dt) - central_freq)
lambda_axis_a = c / (w_axis_a / (2 * np.pi * hbar))

# Experimental Spectra
absorption_exp_130 = np.load(spectra_data + 'absorption_exp_130.npy')

fig, (ax) = plt.subplots(1, 1)

ax.plot(np.fft.fftshift(lambda_axis_a)[9500:10100],
        np.fft.fftshift(Aw_norm)[9500:10100], label="Renger's Hamiltonian",
        linewidth=4.5, color='Blue', alpha=0.6)

ax.plot(absorption_exp_130[:, 0], absorption_exp_130[:, 1] / np.max(
    absorption_exp_130[:, 1]),
        color='slategray', label='Experiment', linewidth=4.5, linestyle='dashed')

ax.set(xlim=(630, 710))
ax.set_xlabel('Wavelength (nm)', fontsize=18)
ax.set_ylabel(f'Absorption (a.u.) at {temp}K, \n N_ens={N_ens},sigma={sigma_e}',
              fontsize=18)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_xticks([630, 640, 660, 680, 700])
ax.set_yticklabels(labels=[0, 0.2, 0.4, 0.6, 0.8, 1], fontsize=15)
ax.set_xticklabels(labels=[630, 640, 660, 680, 700], fontsize=15)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.legend()
plt.tight_layout()
# plt.savefig(f"{working_dir}abso_{temp}K_N_ens_{N_ens}_sigma_{sigma_e}.png")
fig.show()
