import numpy as np
import warnings
from numpy.random import RandomState

def rarefy(x, depth=1000, iterations=1, seed=42):
    """
    Rarefies a count or frequency vector 'x' by randomly subsampling elements.

    Parameters:
    - x (numpy.ndarray): Input count or frequency vector to be rarefied. Meant to represent gene or species counts.
    - depth (int, optional): The desired rarefaction depth, i.e., the number of elements to subsample.
                             Default is 1000.
    - iterations (int, optional): The number of iterations to perform rarefaction.
                                  Default is 1. If > 1, random list of seeds is generated. Overules seed param.
    - seed (int, optional): Seed for reproducibility of random sampling. Default is 42.

    Returns:
    numpy.ndarray: Rarefied vector with the same length as the input vector 'x'.
                  The result is the mean of rarefied counts over multiple iterations.
                  If the number of iterations exceeds 100000, a warning is printed, and the
                  number of iterations is set to 100000.
                  If the rarefaction depth exceeds the total count in 'x', an array of NaNs with the
                  same length as 'x' is returned.
                  If 'x' has zero counts or length zero, an array of NaNs with the same length as 'x' is returned.
    """

    # Convert pandas Series or lists to numpy array if needed
    if not isinstance(x, np.ndarray):
        x = np.array(x)

    noccur = np.sum(x)
    nvar = len(x)

    # Check for invalid count vectors
    if noccur == 0:
        warnings.warn("Input vector has zero counts")
        return np.array([np.nan] * nvar)

    if nvar == 0:
        warnings.warn("Input vector has zero length")
        return np.array([])

    # Check if the number of iterations is within a reasonable range
    if iterations <= 0:
        raise ValueError("Number of iterations must be positive")

    if iterations > 100000:
        warnings.warn(
            'Max number of iterations allowed is 100000, setting to 100000')
        iterations = 100000

    # Check if the rarefaction depth exceeds the total count in 'x'
    if depth <= 0:
        raise ValueError("Rarefaction depth must be positive")

    if depth > noccur:
        warnings.warn(
            f"Rarefaction depth ({depth}) exceeds total count in vector ({noccur})")
        return np.array([np.nan]*nvar)

     # Calculate probability vector
    p = x/noccur
    seeds = np.random.choice(
        100000, size=iterations) if iterations > 1 else [seed]

    # Initialize results array
    results = np.zeros((iterations, nvar))

    # Perform rarefaction for each iteration
    for i, s in enumerate(seeds):
        prng = RandomState(s)
        choice = prng.choice(nvar, size=depth, p=p)
        results[i, :] = np.bincount(choice, minlength=nvar)
    # Return the mean of rarefied counts over multiple iterations
    return np.nanmean(results, axis=0)  # or np.mean?

def rarefy_insertion_file(input_file, output_file, depth=100000, seed=42):
    """
    Rarefy an insertion file to a specific read depth
    
    Parameters:
    - input_file: path to input insertion file (.gz or .tsv)
    - output_file: path to output rarefied file
    - depth: target number of reads
    - seed: random seed for reproducibility
    """
    import pandas as pd
    import gzip
    
    # Load the file
    if input_file.endswith('.gz'):
        with gzip.open(input_file, 'rt') as f:
            df = pd.read_csv(f, sep='\t', header=None, 
                           names=['Chromosome', 'Start', 'End', 'Abundance', 'Strand'])
    else:
        df = pd.read_csv(input_file, sep='\t', header=None,
                       names=['Chromosome', 'Start', 'End', 'Abundance', 'Strand'])
    
    # Rarefy the abundance column
    rarefied_abundance = rarefy(df['Abundance'].values, depth=depth, seed=seed)
    
    # Update the dataframe with rarefied abundances
    df['Abundance'] = rarefied_abundance
    
    # Remove rows with zero abundance after rarefaction
    df = df[df['Abundance'] > 0]
    
    # Round to integers
    df['Abundance'] = df['Abundance'].round().astype(int)
    
    # Save the rarefied file
    df.to_csv(output_file, sep='\t', header=False, index=False)
    
    return df

def rarefy_binomial(x, d, seed=42):
    """
    Rarefy a count vector using binomial sampling.
    
    Parameters:
    - x (numpy.ndarray): Input count vector to be rarefied
    - d (float): Proportion to keep (between 0 and 1)
    - seed (int, optional): Random seed for reproducibility. Default is 42.
    
    Returns:
    numpy.ndarray: Rarefied vector with binomial sampling applied to each element
    """
    # Convert pandas Series or lists to numpy array if needed
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    
    # Validate inputs
    if d < 0 or d > 1:
        raise ValueError("Proportion d must be between 0 and 1")
    
    if len(x) == 0:
        return np.array([])
    
    # Set random seed
    np.random.seed(seed)
    
    # Apply binomial sampling to each element
    # For each count in x, sample from binomial distribution with n=count, p=d
    rarefied = np.array([np.random.binomial(count, d) for count in x])
    
    return rarefied