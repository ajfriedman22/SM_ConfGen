import pandas as pd
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import os

def define_dihedral(top_file : str):
    """
    Parses the topology for the small molecule to determine all potential dihedral angles

    Parameters
    ----------
    Input 
        top : Name of the topolgy file
    Output
        df : Pandas DF defining all unique dihedral angles between heavy atoms in the topology
    """

    #Load topology
    top = open(top_file, 'r').readlines()

    #Determine atom numbers for non-H atoms
    correct_sect = False
    non_H_atoms, atom_name = [], []
    for l, line in enumerate(top):
        if correct_sect:
            line_sep = line.split(' ')
            while '' in line_sep:
                line_sep.remove('')
            if line_sep[0] != ';':
                if line_sep != None and 'H' not in line_sep[4]:
                    non_H_atoms.append(line_sep[0])
                    atom_name.append(line_sep[4])
                if top[l+1] == '\n':
                    correct_sect = False
                    break
        elif 'atoms' in line:
            correct_sect = True
        elif 'bonds' in line:
            correct_sect = False
            break

    #Get unique dihedrals which do not involve hydrogens
    dihe_list = []
    for line in top:
        if correct_sect:
            line_sep = line.split(' ')
            while '' in line_sep:
                line_sep.remove('')
            if line_sep[0] != ';':
                if line_sep != None and line_sep[0] in non_H_atoms and line_sep[1] in non_H_atoms and line_sep[2] in non_H_atoms and line_sep[3] in non_H_atoms and line_sep[0:4] not in dihe_list and line_sep[4:0] not in dihe_list:
                    dihe_list.append(line_sep[0:4])
        elif 'dihedrals' in line:
            correct_sect = True
        elif 'Position' in line:
            correct_sect = False
            break
    dihedral_num, atom_1, atom_2, atom_3, atom_4 = [], [], [], []
    for i, dihe in enumerate(dihe_list):
        dihedral_num.append(i+1)
        atom_1.append(atom_name[int(dihe[0])-1])
        atom_2.append(atom_name[int(dihe[1])-1])
        atom_3.append(atom_name[int(dihe[2])-1])
        atom_4.append(atom_name[int(dihe[3])-1])
    df = pd.DataFrame({'Dihedral Number': dihedral_num, 'Atom Name 1': atom_1, 'Atom Name 2': atom_2, 'Atom Name 3': atom_3, 'Atom Name 4': atom_4})
    return df

def get_dihedral(dihedrals_df : pd.DataFrame, traj : md.trajectory):
    """
    Determines the value for the all dihedral angles in the small molecule

    Parameters
    ----------
    Input
        dihedrals_df : pandas dataframe with dihedral definitions
    Output
        dihedrals : array of dihedral angles
    """
    torsion_ind = np.zeros((len(dihedrals_df.rows), 4))
    torsion_name = []
    for i, row in enumerate(dihedrals_df.rows):
        torsion_name.append(row['Dihedral Number'])
        for j in range(1, 5):
            atom_name = row[f'Atom Name {j}']
            atom = traj.topology.select(f'name {atom_name}')
            if len(atom) > 1:
                raise Exception(f'Error: Multiple atoms with name {atom_name}')
            elif len(atom) == 0:
                raise Exception(f'Error: No atom named {atom_name}')
            torsion_ind[i,j] = atom
    return torsion_name, torsion_ind

def deter_multimodal(dihedrals : np.array, i : int, dihe_name : str):
    """
    Determine the maxima for an n-modal distribution and then plot it

    Parameters

    """
    from scipy.stats import gaussian_kde
    from scipy.signal import find_peaks

    #Seperate dihedral angles
    dihe_dist = dihedrals[:,i]
    
    #Determine modes for distribtion from KDE
    kde_default = gaussian_kde(dihe_dist)
    num_peaks = 0
    for i, x in enumerate(np.logspace(0, kde_default.factor*10, num=100)):
        # Estimate the density function using KDE
        kde = gaussian_kde(dihe_dist, bw_method = x)
        x_vals = np.linspace(min(dihe_dist), max(dihe_dist), 1000)
        density_estimation = kde.evaluate(x_vals)
    
        # Find local maxima in the density function
        peaks, _ = find_peaks(density_estimation)
    
        # Retrieve x values corresponding to the maxima
        maxima = x_vals[peaks]

        # Determine if number of peaks is stable
        if i%10 == 0 and len(maxima) == num_peaks:
            break
        elif i%10 == 0:
            num_peaks = len(maxima)
    
    #Create directory for dihedral distributions if not present
    if not os.path.exists('analysis/dihedrals'):
        os.mkdir('analysis/dihedrals')

    #Histogram of the data
    n, bins, patches = plt.hist(dihe_dist, 100, density=True, facecolor='g', alpha=0.75)
    #Inidcate Maxima
    for i in range(len(maxima)):
        plt.axvline(x = maxima[i], color = 'k')

    plt.xlabel('Torsional Angle(rad)')
    plt.ylabel('Probability')
    plt.xlim(-180, 180)
    plt.title('Histogram of Torsion Angle ' + dihe_name)
    plt.grid(True)
    plt.savefig('analysis/dihedrals/dihe_angle_' + dihe_name + '.png')
    plt.close()
    return maxima

def compute_max(data):
    from scipy.stats import gaussian_kde
    import numpy as np

    kde = gaussian_kde(data)
    samples = np.linspace(min(data), max(data), 50)
    probs = kde.evaluate(samples)
    maxima_index = probs.argmax()
    maxima = samples[maxima_index]
    
    return maxima