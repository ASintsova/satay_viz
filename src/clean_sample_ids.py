def edit_sample_id(sample_id_str):
    """
    First step: Clean and standardize sample ID format.
    
    Args:
        sample_id_str (str): Original sample ID
        
    Returns:
        str: Cleaned sample ID with standardized formatting
    """
    new_id = (
        sample_id_str.replace("-", "_")
        .replace("woaF", "0nMaF")
        .replace("noaF", "0nMaF")
        .replace(".", "_")
    )

    if "20200305" in sample_id_str:
        new_id = new_id.replace("15776_", "3_").replace("20200305_", "")
    else:
        new_id = new_id.replace("20190221_", "")
    
    new_id = new_id.replace("20190223_", "").replace("A_", "rep_")

    return new_id


def edit_mutant_sample_id(sample_id_str):
    """
    Second step: Convert sample ID to final standardized format with mutant names.
    
    Args:
        sample_id_str (str): Sample ID from edit_sample_id function
        
    Returns:
        str: Final formatted sample ID in format "rep_X_mutant_conc"
        
    Raises:
        KeyError: If mutant ID is not found in the dictionary
    """
    mutant_dict = {
        '688': 'whi3dpQ',
        '690': 'ppz1', 
        '693': 'caf40old',
        '16184': 'pin4',
        '16216': '16216',
        '16824': 'whi5',
        '16826': 'ppz1dN',
        '16825': 'cln3',
        '685': 'whi3',
        '693noab': 'caf40',
        '693withab': 'caf40Ab',
        '16825withab': 'cln3Ab'
    }
    
    # Handle rep_ prefix
    if 'rep_' in sample_id_str:
        new_id = sample_id_str.split("rep_")[1]
    else:
        new_id = sample_id_str
    
    # Split into parts
    parts = new_id.split("_")
    
    # Get mutant name
    mutant_key = parts[0]
    if mutant_key not in mutant_dict:
        raise KeyError(f"Unknown mutant ID: {mutant_key}")
    
    mut = mutant_dict[mutant_key]
    
    # Determine replicate and concentration
    if len(parts) == 1:
        # Just the mutant ID
        conc = ''
        rep = 'rep_1'
    elif 'nMaF' in parts[1] or parts[1] in ['0nMaF', '_ab', 'withab']:
        conc = parts[1]
        rep = 'rep_1'
    elif parts[1] == '2' and len(parts) > 2:
        # Format like 693_2_4nMaF -> parts would be ['693', '2', '4nMaF']
        conc = parts[2]
        rep = 'rep_2'
    elif parts[1].endswith('ab') or parts[-1].endswith('ab'):
        rep = 'rep_wAb'
        conc = ''
    else:
        # Default case
        conc = parts[1] if len(parts) > 1 else ''
        rep = 'rep_1'
    
    # Join and remove trailing underscore
    return "_".join([rep, mut, conc]).rstrip("_")


def process_sample_ids(sample_ids):
    """
    Process a list of sample IDs through both functions.
    
    Args:
        sample_ids (list): List of original sample ID strings
        
    Returns:
        list: List of tuples (original, step1_result, final_result)
    """
    results = []
    errors = []
    
    for sample_id in sample_ids:
        try:
            step1_result = edit_sample_id(sample_id)
            final_result = edit_mutant_sample_id(step1_result)
            results.append((sample_id, step1_result, final_result))
        except Exception as e:
            errors.append((sample_id, str(e)))
    
    return results, errors