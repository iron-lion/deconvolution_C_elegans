import numpy as np
import pandas as pd
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

def lsvr_deconvolution(bulk_expression_data, reference_matrix, param_C, param_e):
    """
    Deconvolute bulk RNA-seq data into cell type proportions using 
    Linear Support Vector Regression (SVR), similar to the core algorithm of CIBERSORT.

    This implementation uses a standard SVR fit and then ensures the output 
    proportions are non-negative and sum to 1.

    Args:
        bulk_expression_data (pd.Series or np.array): A single bulk sample's 
                                                     gene expression profile. 
                                                     (Length G, where G is the number of genes)
        reference_matrix (pd.DataFrame or np.array): The signature matrix (X). 
                                                    Rows are genes (G), columns are cell types (C).
                                                    The gene order MUST match bulk_expression_data.

    Returns:
        pd.Series: A series of cell type proportions (C), summing to 1.
    """
    
    # Ensure the reference matrix is a DataFrame for easier handling
    if isinstance(reference_matrix, np.ndarray):
        reference_matrix = pd.DataFrame(reference_matrix)

    # Ensure bulk data is a 1D array/Series
    if isinstance(bulk_expression_data, pd.Series):
        y = bulk_expression_data.values
    else:
        y = np.array(bulk_expression_data)
    
    # Check for matching gene number (rows in reference, elements in bulk)
    if reference_matrix.shape[0] != len(y):
        raise ValueError("The number of genes in the reference matrix rows and the bulk data must match.")

    # X: Reference matrix (G genes x C cell types) -> Features for SVR
    # The SVR requires the input features (X) to be C columns (cell types) and G rows (genes).
    # X_svr will be (G, C)
    X_svr = reference_matrix.values 
    
    # y: Bulk expression (G genes) -> Target variable for SVR
    # y_svr will be (G,)
    y_svr = y

    # Standardize the features (X_svr) and the target (y_svr)
    # Note: Standardizing the features (reference matrix columns) is crucial for SVR performance.
    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X_svr)
    
    # Standardizing the target (y) often improves stability too
    scaler_y = StandardScaler()
    y_scaled = scaler_y.fit_transform(y_svr.reshape(-1, 1)).flatten()
    

    svr = SVR(kernel='linear', C=param_C, epsilon=param_e)
    svr.fit(X_scaled, y_scaled)
    
    # The coefficients (weights) are the cell type proportions
    proportions = svr.coef_[0]

    # A. Enforce Non-Negativity
    # Set any negative proportion to zero, as cell type abundance cannot be negative.
    proportions = proportions.copy()
    proportions[proportions < 0] = 0
    
    # B. Normalize to Sum to 1
    # If the sum is already 0 (e.g., all were negative), avoid division by zero
    total_proportion = np.sum(proportions)
    if total_proportion > 0:
        normalized_proportions = proportions / total_proportion
    else:
        normalized_proportions = np.zeros_like(proportions)
        
    cell_type_names = reference_matrix.columns if isinstance(reference_matrix, pd.DataFrame) else [f'CellType_{i+1}' for i in range(proportions.shape[0])]
    
    result_series = pd.Series(normalized_proportions, index=cell_type_names)
    
    return result_series


def run_multiprocess_deconvolution(bulk_matrix, reference_matrix, param_C=0.001, param_e=0.1, max_workers=None):
    """
    Parallelizes the SVR deconvolution across all bulk samples.

    Args:
        bulk_matrix (pd.DataFrame): DataFrame of bulk samples (Samples x Genes).
        reference_matrix (pd.DataFrame): Signature matrix (Genes x Cell Types).
        max_workers (int, optional): The maximum number of processes to use. 
                                     Defaults to the number of CPU cores.

    Returns:
        list[pd.Series]: A list of cell type proportions for each bulk sample.
    """
    # Create an iterable of individual bulk samples (rows) to feed to the map function.
    # Note the .T (transpose) is applied to match the expected format (Genes as index).
    bulk_samples_iterable = [bulk_matrix.iloc[i,:].T for i in range(len(bulk_matrix))]

    # Use functools.partial to fix the 'reference_matrix' argument for the worker function.
    # This is crucial because the map function only iterates over the first argument.
    deconvolute_partial = partial(lsvr_deconvolution, reference_matrix=reference_matrix, param_C=param_C, param_e=param_e)

    # Use ProcessPoolExecutor to manage the worker processes
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # executor.map applies the function (deconvolute_partial) to every item 
        # in the iterable (bulk_samples_iterable) in parallel.
        
        final_proportions = list(tqdm(
            executor.map(deconvolute_partial, bulk_samples_iterable), 
            total=len(bulk_samples_iterable), 
            desc="Deconvoluting Samples"
        ))

    return final_proportions
