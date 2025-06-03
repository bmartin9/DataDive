import numpy as np

def pg_rapid_labf(Zsq, Wsq):
    return 0.5 * (np.log(Wsq) + Zsq * (1 - Wsq))

def pg_get_labfs(effect_est, effect_se, binary_outcomes):
    effect_est = np.asarray(effect_est)
    effect_se = np.asarray(effect_se)
    binary_outcomes = np.asarray(binary_outcomes)

    if np.isnan(effect_est).any():
        raise ValueError("there are missing values in effect.est")
    if np.isnan(effect_se).any():
        raise ValueError("there are missing values in effect.se")
    if np.any(effect_se == 0):
        raise ValueError("there are zero values in effect.se")
    if np.isnan(binary_outcomes).any():
        raise ValueError("there are missing values in binary.outcomes")
    if not np.all(np.isin(binary_outcomes, [0, 1])):
        raise ValueError("binary.outcomes should contain values that are either 0s or 1s")

    if effect_est.ndim == 2:
        if len(binary_outcomes) != effect_est.shape[1]:
            raise ValueError("binary.outcomes should match the number of columns in effect.est")
    else:
        if binary_outcomes.size != 1:
            raise ValueError("binary.outcomes should be a single value or a vector of length 1")

    Z = effect_est / effect_se
    W = 0.15 / effect_se

    if effect_se.ndim == 2:
        for idx, is_binary in enumerate(binary_outcomes):
            if is_binary == 1:
                W[:, idx] = 0.2 / effect_se[:, idx]
    else:
        if binary_outcomes == 1:
            W = 0.2 / effect_se

    Zsq = Z ** 2
    Wsq = 1 / (1 + W ** 2)

    return pg_rapid_labf(Zsq, Wsq)

if __name__ == "__main__":
    import pandas as pd

    DATA_DIR = 'C:/Users/alexy/Downloads/simulated_data_bioxcelerate/DataDive/data/processed/concatenated_data_csv/split_0/'
    beta_file_names = [
        'betas_train.csv',
        'betas_validate.csv',
        'betas_test.csv']
    se_file_names = [
        'se_train.csv',
        'se_validate.csv',
        'se_test.csv']
    destinations = [
        'labfs_train.csv',
        'labfs_validate.csv',
        'labfs_test.csv']
    for i in range(3):
        betas = pd.read_csv(DATA_DIR + beta_file_names[i])
        se = pd.read_csv(DATA_DIR + se_file_names[i])

        # Ensure the 'file_path' column is present
        if 'file_path' not in betas.columns or 'file_path' not in se.columns:
            raise ValueError("The 'file_path' column is missing from the input data.")
        betas_vals = betas[[col for col in betas.columns if col.startswith('beta')]].values
        se_vals = se[[col for col in se.columns if col.startswith('se')]].values
        binary_outcomes = np.zeros(betas_vals.shape[1], dtype=int)  # Assuming binary outcomes are all zeros for this example
        labfs = pg_get_labfs(betas_vals, se_vals, binary_outcomes)

        # Create labf column names
        labf_cols = [f'labf_{i}' for i in range(labfs.shape[1])]

        # Build the DataFrame
        labfs_df = pd.DataFrame(
            labfs,
            columns=labf_cols
        )
        labfs_df.insert(0, 'file_path', betas['file_path'].values)
        labfs_df['cv_idx'] = betas['cv_idx'].values
        labfs_df['cluster_num'] = betas['cluster_num'].values

        # Reorder columns: file_path, labf_0...labf_n, cv_idx, cluster_num
        cols = ['file_path'] + labf_cols + ['cv_idx', 'cluster_num']
        labfs_df = labfs_df[cols]

        # Save to CSV if needed
        labfs_df.to_csv(DATA_DIR + destinations[i], index=False)


