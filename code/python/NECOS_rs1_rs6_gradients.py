import argparse
import os
import sys

# Constants
NOcomp = 10
kernel = 'normalized_angle'
approach = 'pca'
INPUT_DATA_ROOT = '/path/to/input_data/rsFC'
OUTPUT_ROOT = '/path/to/output_data/multi_rest_gradients'

parser = argparse.ArgumentParser(description='NECOS resting-state gradient analysis')
parser.add_argument(
    '--dry-run',
    action='store_true',
    help='Validate input and output paths without processing any files.',
)
args = parser.parse_args()

# Create subject list for all 67 participants
sub = [f'NECOS{str(i).zfill(3)}' for i in range(1, 68)]
resting_states = range(1,7)

# Main directory
main_dir = OUTPUT_ROOT

resting_states = [6]

if args.dry_run:
    print(f'Input data root: {INPUT_DATA_ROOT}')
    print(f'Output root: {OUTPUT_ROOT}')

    has_error = False
    if INPUT_DATA_ROOT.startswith('/path/to/'):
        print('INPUT_DATA_ROOT still uses a placeholder path.')
        has_error = True
    if OUTPUT_ROOT.startswith('/path/to/'):
        print('OUTPUT_ROOT still uses a placeholder path.')
        has_error = True

    if not os.path.exists(INPUT_DATA_ROOT):
        print('Input data root does not exist.')
        has_error = True
    else:
        for rs in resting_states:
            available_files = 0
            for s in sub:
                fmri_filename = os.path.join(
                    INPUT_DATA_ROOT,
                    s,
                    'functional',
                    f'rest{rs}',
                    'rest_preprocessed_hp_filt_2mni_6fwhm.nii.gz',
                )
                if os.path.exists(fmri_filename):
                    available_files += 1
            print(f'rest{rs}: found {available_files}/{len(sub)} expected files')

    if os.path.exists(OUTPUT_ROOT):
        print('Output root already exists and will be reused.')
    else:
        print('Output root does not exist yet and will be created during a normal run.')

    sys.exit(1 if has_error else 0)

import numpy as np
import pandas as pd
from nilearn.datasets import fetch_atlas_schaefer_2018
from nilearn.maskers import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from brainspace.gradient import GradientMaps
from brainspace.datasets import load_group_fc


# Fetch atlas and create masker
dataset = fetch_atlas_schaefer_2018(n_rois=400, resolution_mm=2)
atlas_filename = dataset.maps
masker = NiftiLabelsMasker(labels_img=atlas_filename, standardize=True)


# Availability DataFrame
availability_df = pd.DataFrame(index=sub, columns=[f'rest{rs}' for rs in range(1, 7)])

# Process each resting state and subject
for rs in resting_states:
    rs_folder = os.path.join(main_dir, f'rs_{rs}')
    os.makedirs(rs_folder, exist_ok=True)

    group_ts = []
    available_sub = []  # Keep track of available subjects for this resting state

    for s in sub:
        fmri_filename = os.path.join(
            INPUT_DATA_ROOT,
            s,
            'functional',
            f'rest{rs}',
            'rest_preprocessed_hp_filt_2mni_6fwhm.nii.gz',
        )

        # Check if file exists
        if os.path.exists(fmri_filename):
            availability_df.at[s, f'rest{rs}'] = 'Yes'
            ts = masker.fit_transform(fmri_filename)
            group_ts.append(ts)
            available_sub.append(s)  # Add subject to available list

            # Save time series in participant-specific folder
            participant_folder = os.path.join(rs_folder, s)
            os.makedirs(participant_folder, exist_ok=True)
            np.save(os.path.join(participant_folder, f'timeseries_hp_filt_2mni_6fwhm_timeseries.npy'), ts)
        else:
            availability_df.at[s, f'rest{rs}'] = 'No'

    # If there are time series data for this resting state
    if group_ts:
        # Connectivity matrices for the group
        correlation_measure = ConnectivityMeasure(kind='correlation')
        conn_matrix_ind = correlation_measure.fit_transform([group_ts][0])
        conn_matrix_ind_fisher = np.arctanh([conn_matrix_ind][0])
        conn_matrix_ind_fisher[np.isposinf([conn_matrix_ind_fisher][0])] = 1

        # Calculate and save group connectivity and affinity matrix
        conn_matrix_group_fisher = np.mean(conn_matrix_ind_fisher, axis=0)
        conn_matrix_group_fisher[np.isposinf(conn_matrix_group_fisher)] = 1

        # load HCP data from BrainSpace and calculate group level gradients
        conn_matrix = load_group_fc('schaefer', scale=400)
        ref_gradient = GradientMaps(n_components=NOcomp, kernel=kernel, approach=approach)
        ref_gradient.fit(conn_matrix)

        group_gradient = GradientMaps(n_components=NOcomp, kernel=kernel, approach=approach, alignment='procrustes')
        group_gradient.fit([conn_matrix_group_fisher], reference=ref_gradient.gradients_)

        np.savetxt(os.path.join(rs_folder, 'group_grad_pca_hcp_aligned.csv'), group_gradient.aligned_[0], delimiter=",")

        # Calculate and save individual gradients
        for i, s in enumerate(available_sub):
            ind_gradient = GradientMaps(n_components=NOcomp, kernel=kernel, approach=approach, alignment='procrustes')
            ind_gradient.fit(conn_matrix_ind_fisher[i], reference=group_gradient.aligned_[0])

            participant_folder = os.path.join(rs_folder, s)
            np.savetxt(os.path.join(participant_folder, f'ind_grad_pca_{s}.csv'), ind_gradient.aligned_, delimiter=",")


# Save availability data
availability_df.to_csv(os.path.join(main_dir, 'data_availability.csv'))