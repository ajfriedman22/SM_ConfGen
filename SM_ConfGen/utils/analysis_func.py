import pandas as pd
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import os
import seaborn as sns

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
    dihedral_num, dihedral_atom_indx = [], []
    for i, dihe in enumerate(dihe_list):
        dihedral_num.append(i+1)
        dihedral_atom_indx.append([int(dihe[0])-1, int(dihe[1])-1, int(dihe[2])-1, int(dihe[3])-1])
    return dihedral_num, dihedral_atom_indx

def get_dihedral(dihedrals_df : pd.DataFrame, traj : md.Trajectory):
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
    
    #Determine maxima for probability distribution
    maxima, prob_minima, width = compute_max(dihe_dist, True)
    
    #If all angles are sampled this much than this is not a multimodal peak
    if prob_minima > 0.002:
        return [maxima], dihe_dist

    #Determine data not in the main peak
    main_peak, other_peak = [], []
    for i in dihe_dist:
        if abs(i - maxima) < width or abs(i + 360 - maxima) < width or abs(i - 360 - maxima) < width:
            main_peak.append(i)
        else:
            other_peak.append(i)
    all_maxima = [maxima]

    #If greater than 5% outliers count as seperate peak
    while len(other_peak)/len(dihe_dist) > 0.05:
        maxima, width = compute_max(other_peak)
        new_dist = other_peak
        main_peak, other_peak = [], []
        for i in new_dist:
            if abs(i - maxima) < width or abs(i + 360 - maxima) < width or abs(i - 360 - maxima) < width:
                main_peak.append(i)
            else:
                other_peak.append(i)
        if len(main_peak) > (0.10*len(dihe_dist)):
            all_maxima.append(maxima)
        else:
            break
    return all_maxima, dihe_dist

def plot_torsion(dihe_dist, maxima, dihe_name):
    #Histogram of the data
    n, bins, patches = plt.hist(dihe_dist, 100, density=True, facecolor='g', alpha=0.75)
    #Inidcate Maxima
    for i in range(len(maxima)):
        plt.axvline(x = maxima[i], color = 'k')

    plt.xlabel('Torsional Angle(rad)')
    plt.ylabel('Probability')
    plt.xlim(-180, 180)
    plt.title(f'Histogram of Torsion Angle {dihe_name}')
    plt.grid(True)
    plt.savefig(f'analysis/dihedrals/dihe_angle_{dihe_name}.png')
    plt.close()

def compute_max(data, init=False):
    from scipy.stats import gaussian_kde
    import numpy as np

    kde = gaussian_kde(data)
    samples = np.linspace(-180, 180, 360)
    probs = kde.evaluate(samples)
    maxima_index = probs.argmax()
    prob_minima = probs.min()
    maxima = samples[maxima_index]

    #Determine peak width
    prob = probs[maxima_index]
    width = 1
    index = maxima_index
    while prob > 0.0005:
        width += 1
        if maxima_index < 100:
            index += 1
        else:
            index -= 1
        if index == 360 or index < 0:
            break
        prob = probs[index]
    if init:
        return maxima, prob_minima, width
    else:
        return maxima, width

def input_torsion():
    df_max = pd.read_csv('analysis/dihe_ind_max.csv')
    torsion_name = list(set(df_max['Dihedral Name'].to_list()))
    torsion_ind = np.zeros((len(torsion_name), 4))
    max_values, peak_options = [], []
    for i, torsion in enumerate(torsion_name):
        df = df_max[df_max['Dihedral Name'] == torsion]
        torsion_ind[i, :] = [int(df.iloc[0, 2]), int(df.iloc[0, 3]), int(df.iloc[0, 4]), int(df.iloc[0, 5])]
        maxima = df['Maxima'].to_list()
        max_values.append(maxima)
        opt = np.linspace(1, len(maxima), num=len(maxima), dtype=int)
        peak_options.append(list(opt))
    return torsion_name, torsion_ind, max_values, peak_options

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    idx2 = (np.abs(array - value + 360)).argmin()
    idx3 = (np.abs(array - value - 360)).argmin()
    if np.abs(array[idx]-value) <= np.abs(array[idx2] - value + 360) and np.abs(array[idx] - value) <= np.abs(array[idx3] - value - 360):
        return idx+1
    elif np.abs(array[idx2] - value + 360) <= np.abs(array[idx] - value) and np.abs(array[idx2] - value + 360) <= np.abs(array[idx3] - value - 360):
        return idx2+1
    else:
        return idx3+1

def clust_conf(traj, per, file_name):
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
    
    #Compute pairwise RMSD
    distances = np.empty((traj.n_frames, traj.n_frames))
    for i in range(traj.n_frames):
        distances[i] = md.rmsd(traj, traj, i, atom_indices=traj.topology.select('element != H'))
    
    #Perform Clustering
    reduced_distances = squareform(distances, checks=False)
    #np.savetxt('test_reduce_dist.txt', reduced_distances)
    link = linkage(reduced_distances, method='single') #The hierarchical clustering encoded as a matrix
    frame_list = dendrogram(link, no_labels=False, count_sort='descendent')['leaves']
    frame_cat = dendrogram(link, no_labels=False, count_sort='descendent')['color_list']

    #Keep only one file per cluster
    frames_sep = [] #List of frames that are unique and will be processed
    cat = frame_cat[0]
    frames_indv = [frame_list[0]]
    for frame in range(1, len(frame_list)-1):
        if frame_cat[frame] == cat:
            frames_indv.append(frame_list[frame])
        else:
            frames_sep.append(frames_indv)
            cat = frame_cat[frame]
            frames_indv = [frame_list[frame]]
    frames_sep.append(frames_indv)
    
    #Analyze each cluster (determine centroid and plot Intercluster RMSD)
    frames_unique, per_unique= compare_within_cluster(traj, frames_sep, per, distances[0], file_name)

    return frames_unique, per_unique, frames_sep

def compare_within_cluster(traj, cluster_frames, per, rmsd_all, file_name):
    per_unique = np.zeros(len(cluster_frames))
    frames_unique = []
    df_clust = pd.DataFrame()
    for f, frames in enumerate(cluster_frames):
        cluster_traj = traj.slice(frames)
        rmsd_clust = np.empty((cluster_traj.n_frames, cluster_traj.n_frames))
        rmsd_clust_mean = np.zeros(cluster_traj.n_frames)
        for i in range(cluster_traj.n_frames):
            rmsd_clust[i] = md.rmsd(cluster_traj, cluster_traj, i, atom_indices=cluster_traj.topology.select('element != H'))
            rmsd_clust_mean[i] = np.mean(rmsd_clust[i])
        min_index = np.argmin(rmsd_clust_mean)
        centroid_frame = frames[min_index]
        frames_unique.append(centroid_frame)
        for frame in frames:
            per_unique[f] += per[frame]
        df = pd.DataFrame({'Cluster ID': f+1, r'RMSD ($\AA$)': rmsd_clust[min_index]*10})
        df_clust = pd.concat([df_clust, df])
    #Save CSV
    df_clust.to_csv(f'analysis/{file_name}-clust-rmsd.csv')

    #Plot Comparison
    plt.figure()
    g = sns.FacetGrid(df_clust, col="Cluster ID", col_wrap=5, xlim=(0,3))
    g.map(sns.histplot, r'RMSD ($\AA$)')
    plt.savefig(f'analysis/{file_name}-clust-rmsd.png')
    plt.close()    
    return frames_unique, per_unique

def process_confs(raw_traj, frames, per, file_name, id_type='Conformer'):
    per_non_zero = np.array(per[per!=0], dtype=float)
    ordered_index = per_non_zero.argsort()
    per_ordered, frames_ordered = [], []
    for idx in reversed(ordered_index):
        per_ordered.append(per[idx])
        frames_ordered.append(frames[idx])
    #Save PDB
    traj = raw_traj.slice(frames_ordered)
    traj.save_pdb(f'{file_name}.pdb')

    #Compute relative conformer energy
    rel_ener = get_rel_ener(per_ordered)

    #compute radius of gyration
    rg = md.compute_rg(traj)

    #Save CSV
    df_clust = pd.DataFrame({f'{id_type} ID': np.linspace(1, len(per), num=len(per), dtype=int), 'Occupancy': per_ordered, 'Relative FE': rel_ener})
    df_clust['Radius of Gyration'] = rg
    df_clust.to_csv(f'analysis/{file_name}.csv')

    labels = []
    for i, per in enumerate(df_clust['Occupancy']):
        if per > 1.5:
            labels.append(df_clust[f'{id_type} ID'].values[i])
        else:
            labels.append('')
 
    plt.figure()
    plt.pie(df_clust['Occupancy'], labels=labels)
    plt.savefig(f'analysis/{file_name}_pie.png')
    plt.close()

    plt.figure()
    sns.barplot(df_clust, x=f'{id_type} ID', hue='Relative FE', y='Relative FE', palette='cool_r', legend=False, order=df_clust[f'{id_type} ID'])
    plt.xlabel('Cluster ID', fontsize=16)
    plt.ylabel('Relative Free Energy (kcal/mol)', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)
    plt.savefig(f'analysis/{file_name}_FE.png')
    plt.close()

    plt.figure()
    sns.barplot(df_clust, x=f'{id_type} ID', y='Radius of Gyration')
    plt.ylabel(r'Radius of Gyration($\AA$)')
    plt.savefig(f'analysis/{file_name}_rg.png')
    plt.close()

    plt.figure()
    sns.histplot(df_clust, x='Radius of Gyration')
    plt.xlabel(r'Radius of Gyration($\AA$)')
    plt.savefig(f'analysis/{file_name}_rg_hist.png') 
    plt.close()

    return ordered_index+1, per.non_zero()[0], traj

def get_rel_ener(per_all):
    from scipy.constants import k, Avogadro
    ref_per = np.max(per_all)
    rel_ener = np.zeros(len(per_all))
    for p, per in enumerate(per_all):
        rel_ener[p] = -(k/1000)*300*np.log(per/ref_per)*Avogadro/4.184
    return rel_ener

