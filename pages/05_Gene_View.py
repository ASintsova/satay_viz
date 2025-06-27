import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import gzip
import numpy as np
import tempfile
import os
import sys
sys.path.append('src')
from clean_sample_ids import edit_sample_id

def clean_filename(filename):
    """Clean up filename using only the edit_sample_id function"""
    try:
        # Extract base name and split at .bam
        base_name = filename.split(".bam")[0]
        # Apply only the first cleanup step
        cleaned_name = edit_sample_id(base_name)
        return cleaned_name
    except Exception as e:
        # If cleanup fails, return original filename
        st.warning(f"Could not clean filename '{filename}': {e}")
        return filename

@st.cache_data
def get_annotations(annotation_file="annotations.tsv"):
    """Load gene annotations with caching"""
    return pd.read_table(annotation_file)

@st.cache_data
def get_chromosome_mapping(chr_file="chr_sizes.tsv"):
    """Load chromosome mapping between chr names and RefSeq IDs with caching"""
    chr_map = pd.read_table(chr_file)
    # Create mapping dict: chrXII -> NC_001144.5
    mapping = dict(zip(chr_map['Chromosome'], chr_map['RefSeq ID']))
    return mapping

def load_insertion_file(file_path):
    """Load insertion file (gzipped or not) into DataFrame"""
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', header=None, 
                               names=['Chromosome', 'Start', 'End', 'Abundance', 'Strand'])
        else:
            df = pd.read_csv(file_path, sep='\t', header=None,
                           names=['Chromosome', 'Start', 'End', 'Abundance', 'Strand'])
        return df
    except Exception as e:
        st.error(f"Error loading file {file_path}: {e}")
        return None

def find_overlapping_insertions_batch(gene_info, file_data_dict, chr_mapping, padding=500):
    """Find overlapping insertions for all files at once - much faster"""
    if gene_info.empty:
        return {}
    
    gene_chr = gene_info['Chromosome'].iloc[0]
    gene_start = gene_info['Start'].iloc[0] - padding
    gene_end = gene_info['End'].iloc[0] + padding
    
    # Map chromosome name to RefSeq ID
    refseq_chr = chr_mapping.get(gene_chr)
    if refseq_chr is None:
        st.warning(f"No RefSeq mapping found for chromosome {gene_chr}")
        return {}
    
    overlapping_results = {}
    
    # Process all files in batch
    for file_name, insertions_df in file_data_dict.items():
        if insertions_df.empty:
            overlapping_results[file_name] = pd.DataFrame()
            continue
            
        # Vectorized filtering - much faster than individual operations
        mask = (
            (insertions_df['Chromosome'] == refseq_chr) &
            (insertions_df['Start'] <= gene_end) & 
            (insertions_df['End'] >= gene_start)
        )
        
        overlapping_results[file_name] = insertions_df[mask].copy()
    
    return overlapping_results


def create_combined_gene_visualization(gene_info, overlapping_data_dict, padding=500):
    """Create combined visualization for all files at once - much faster"""
    if gene_info.empty:
        return go.Figure(), []
    
    gene_start = gene_info['Start'].iloc[0]
    gene_end = gene_info['End'].iloc[0]
    gene_name = gene_info['gene'].iloc[0]
    gene_strand = gene_info['Strand'].iloc[0]
    
    # Expanded region coordinates
    region_start = gene_start - padding
    region_end = gene_end + padding
    
    fig = go.Figure()
    file_stats = []
    
    # Calculate global abundance range for consistent coloring
    all_abundances = []
    for overlapping_insertions in overlapping_data_dict.values():
        if not overlapping_insertions.empty:
            all_abundances.extend(overlapping_insertions['Abundance'].values)
    
    if all_abundances:
        global_min_log = np.log10(min(all_abundances) + 1)
        global_max_log = np.log10(max(all_abundances) + 1)
    else:
        global_min_log = global_max_log = 0
    
    # Sort files for consistent display
    sorted_files = sorted(overlapping_data_dict.keys())
    
    # Process all files at once
    for i, file_name in enumerate(sorted_files):
        overlapping_insertions = overlapping_data_dict[file_name]
        y_offset = i * 2
        
        # Clean filename for display
        display_name = clean_filename(file_name)
        
        # Draw expanded region background
        fig.add_shape(
            type="rect",
            x0=region_start, y0=y_offset-0.8, x1=region_end, y1=y_offset+0.8,
            line=dict(color="lightblue", width=1),
            fillcolor="lightblue",
            opacity=0.2
        )
        
        # Draw gene rectangle (transparent fill so insertion lines show through)
        fig.add_shape(
            type="rect",
            x0=gene_start, y0=y_offset-0.5, x1=gene_end, y1=y_offset+0.5,
            line=dict(color="black", width=2),
            fillcolor="rgba(255,255,255,0.3)",  # Semi-transparent white
            opacity=0.8,
            layer="below"  # Ensure rectangle is below traces
        )
        
        # Add file name annotation
        fig.add_annotation(
            x=region_start,
            y=y_offset+0.6,
            text=display_name,
            showarrow=False,
            font=dict(size=8, color="darkblue"),
            xanchor="left"
        )
        
        # Add strand direction indicator
        if gene_strand == '+':
            fig.add_annotation(
                x=region_start, y=y_offset-0.6, text="5'",
                showarrow=False, font=dict(size=10, color="black"), xanchor="left"
            )
        elif gene_strand == '-':
            fig.add_annotation(
                x=region_end, y=y_offset-0.6, text="5'",
                showarrow=False, font=dict(size=10, color="black"), xanchor="right"
            )
        
        # Collect file stats
        total_insertions = len(overlapping_data_dict[file_name]) if file_name in overlapping_data_dict else 0
        stats = {
            'File': display_name,
            'Overlapping Insertions': len(overlapping_insertions)
        }
        file_stats.append(stats)
        
        # Add insertion lines using vectorized scatter plot - much faster than individual shapes
        if not overlapping_insertions.empty:
            # Vectorized calculations
            x_positions = (overlapping_insertions['Start'] + overlapping_insertions['End']) / 2
            log_abundances = np.log10(overlapping_insertions['Abundance'] + 1)
            
            # Batch color calculations
            if global_max_log > global_min_log:
                normalized_abundances = (log_abundances - global_min_log) / (global_max_log - global_min_log)
            else:
                normalized_abundances = np.full(len(log_abundances), 0.5)
            
            # Batch approach: Group insertions by color intensity to reduce number of traces
            # Bin normalized abundances into 10 color groups for efficiency
            n_color_bins = 10
            color_bins = np.linspace(0, 1, n_color_bins + 1)
            
            for bin_idx in range(n_color_bins):
                bin_min = color_bins[bin_idx]
                bin_max = color_bins[bin_idx + 1]
                
                # Find insertions in this color bin
                mask = (normalized_abundances >= bin_min) & (normalized_abundances < bin_max)
                if bin_idx == n_color_bins - 1:  # Include the maximum value in the last bin
                    mask = (normalized_abundances >= bin_min) & (normalized_abundances <= bin_max)
                
                if not mask.any():
                    continue
                
                bin_x_positions = x_positions[mask]
                bin_abundances = overlapping_insertions['Abundance'][mask]
                bin_log_vals = log_abundances[mask]
                
                # Calculate color for this bin
                bin_color = (bin_min + bin_max) / 2
                red_val = int(bin_color * 255)
                blue_val = int((1 - bin_color) * 255)
                color = f"rgb({red_val},0,{blue_val})"
                
                # Create line data for this color group
                line_x = []
                line_y = []
                hover_texts = []
                
                for x_pos, abundance, log_val in zip(bin_x_positions, bin_abundances, bin_log_vals):
                    # Add line coordinates: fit within gene rectangle
                    line_x.extend([x_pos, x_pos, None])
                    line_y.extend([y_offset-0.4, y_offset+0.4, None])  # Fit within gene rect (±0.5)
                    hover_texts.extend([
                        f"File: {display_name}<br>Abundance: {abundance}<br>Log: {log_val:.2f}",
                        "", ""  # Only first point has hover
                    ])
                
                # Add this color group as a single trace
                fig.add_trace(go.Scatter(
                    x=line_x,
                    y=line_y,
                    mode='lines',
                    line=dict(width=3, color=color),  # Increased width for visibility
                    hovertemplate='%{text}<extra></extra>',
                    text=hover_texts,
                    showlegend=False,
                    name=f"insertions_{i}_bin_{bin_idx}",
                    opacity=0.9  # Slightly transparent to see overlapping lines
                ))
    
    # Set up layout
    num_files = len(sorted_files)
    y_min = -0.8 - 0.5
    y_max = (num_files-1)*2 + 0.8 + 0.5
    
    fig.update_layout(
        title=f"Gene View: {gene_name} (±{padding}bp)",
        xaxis_title="Genomic Position",
        yaxis_title="Files",
        yaxis=dict(range=[y_min, y_max], showticklabels=False),
        height=200 + num_files * 150,
        showlegend=False,
        template="plotly_white",
        xaxis=dict(range=[region_start, region_end])
    )
    
    # Add colorbar if we have data
    if global_min_log < global_max_log:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                colorscale=[[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']],
                cmin=global_min_log,
                cmax=global_max_log,
                colorbar=dict(title="Log10(Abundance + 1)")
            ),
            showlegend=False
        ))
    
    return fig, file_stats

def get_original_read_count(file_path_or_uploaded_file, is_uploaded=False):
    """Helper function to get original read count from a file"""
    try:
        if is_uploaded:
            # Handle uploaded file
            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_path_or_uploaded_file.name.split('.')[-1]}") as tmp_file:
                tmp_file.write(file_path_or_uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            df = load_insertion_file(tmp_file_path)
            os.unlink(tmp_file_path)
        else:
            # Handle file path
            df = load_insertion_file(file_path_or_uploaded_file)
        
        return int(df['Abundance'].sum()) if df is not None else 0
    except Exception:
        return 0

def apply_rarefaction_to_dataframe(df, depth, total_reads, seed=42):
    """Apply binomial rarefaction to a dataframe and return processed dataframe"""
    from rarefy import rarefy_binomial
    
    if df.empty:
        return df
    
    try:
        # Calculate proportion to keep
        d = min(1.0, depth / total_reads)  # Cap at 1.0 if depth > total_reads
        
        # Apply binomial rarefaction
        rarefied_abundance = rarefy_binomial(df['Abundance'].values, d=d, seed=seed)
        df = df.copy()
        df['Abundance'] = rarefied_abundance
        
        # Remove zero abundance entries
        df = df[df['Abundance'] > 0]
        
        return df
    except Exception as e:
        st.error(f"Error during rarefaction: {e}")
        return df

def get_cache_key(file_names, depth, use_rarefaction):
    """Generate cache key for processed files"""
    file_key = "_".join(sorted(file_names))
    return f"{file_key}_{depth if use_rarefaction else 'no_rarefy'}"


def process_files_fresh(uploaded_files, use_rarefaction, rarefaction_depth, original_read_counts):
    """Process files without caching - actual file processing logic"""
    file_data = {}
    
    if not uploaded_files:
        # Use example files
        example_files = [
            "examples/20190221.A-2_noaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz",
            "examples/20190221.A-2_4nMaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz"
        ]
        
        for file_path in example_files:
            file_name = file_path.split('/')[-1]
            insertions_df = load_insertion_file(file_path)
            
            if insertions_df is not None:
                total_reads = original_read_counts[file_name]
                
                if use_rarefaction:
                    insertions_df = apply_rarefaction_to_dataframe(insertions_df, rarefaction_depth, total_reads)
                    processed_name = f"{file_name}_rarefied_{rarefaction_depth}"
                else:
                    processed_name = file_name
                
                file_data[processed_name] = insertions_df
    else:
        # Process uploaded files
        for uploaded_file in uploaded_files:
            # Create temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{uploaded_file.name.split('.')[-1]}") as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            
            insertions_df = load_insertion_file(tmp_file_path)
            os.unlink(tmp_file_path)
            
            if insertions_df is not None:
                total_reads = original_read_counts[uploaded_file.name]
                
                if use_rarefaction:
                    insertions_df = apply_rarefaction_to_dataframe(insertions_df, rarefaction_depth, total_reads)
                    processed_name = f"{uploaded_file.name}_rarefied_{rarefaction_depth}"
                else:
                    processed_name = uploaded_file.name
                
                file_data[processed_name] = insertions_df
    
    return file_data

def calculate_read_counts_and_depths(uploaded_files):
    """Calculate total read counts for all files and determine max/default depths"""
    read_counts = {}
    
    if not uploaded_files:
        # Use example files
        example_files = [
            "examples/20190221.A-2_noaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz",
            "examples/20190221.A-2_4nMaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz"
        ]
        for file_path in example_files:
            file_name = file_path.split('/')[-1]
            total_reads = get_original_read_count(file_path, is_uploaded=False)
            read_counts[file_name] = total_reads
    else:
        # Calculate for uploaded files
        for uploaded_file in uploaded_files:
            total_reads = get_original_read_count(uploaded_file, is_uploaded=True)
            read_counts[uploaded_file.name] = total_reads
    
    if not read_counts:
        return {}, 1000000, 100000
    
    # Calculate max and default depths
    min_reads = min(read_counts.values())
    max_depth = min_reads  # Max depth is the smallest file's total reads
    default_depth = max(1000, min_reads // 10)  # Default is 1/10 of smallest, min 1000
    
    return read_counts, max_depth, default_depth

def app():
    st.header("Gene-Level Insertion Visualization")
    
    # Load annotations and chromosome mapping
    annotations = get_annotations()
    chr_mapping = get_chromosome_mapping()
    
    # Gene selection
    st.sidebar.header("Gene Selection")
    available_genes = sorted(annotations['gene'].dropna().unique())
    
    # Default to STE11 if available
    default_gene = "STE11" if "STE11" in available_genes else available_genes[0]
    selected_gene = st.sidebar.selectbox("Select gene:", available_genes, 
                                       index=available_genes.index(default_gene))
    
    # File upload and rarefaction controls
    st.sidebar.header("Data Upload")
    uploaded_files = st.sidebar.file_uploader("Upload insertion files (.gz or .tsv)", 
                                            type=['gz', 'tsv', 'txt'], 
                                            accept_multiple_files=True)
    
    # Calculate read counts and determine depth limits
    original_read_counts, max_rarefaction, default_rarefaction = calculate_read_counts_and_depths(uploaded_files)
    
    # Display file information
    if uploaded_files:
        st.sidebar.subheader("Original File Information")
        for uploaded_file in uploaded_files:
            total_reads = original_read_counts[uploaded_file.name]
            st.sidebar.text(f"{uploaded_file.name}: {total_reads:,} reads")
    else:
        st.sidebar.subheader("Example File Information")
        for file_name, total_reads in original_read_counts.items():
            st.sidebar.text(f"{file_name}: {total_reads:,} reads")
    
    st.sidebar.header("Rarefaction Settings")
    use_rarefaction = st.sidebar.checkbox("Apply rarefaction", value=True)
    
    rarefaction_depth = st.sidebar.number_input(
        f"Rarefaction depth (max: {max_rarefaction:,})", 
        min_value=1000, 
        max_value=max_rarefaction, 
        value=default_rarefaction, 
        step=1000
    )
    
    st.sidebar.info(f"Using binomial rarefaction: each file downsampled proportionally to reach target depth")
    
    # Load and process files with caching (only when rarefaction settings change)
    # Cache key based on files and rarefaction settings, NOT gene selection
    if uploaded_files:
        file_names_for_cache = [f.name for f in uploaded_files]
    else:
        file_names_for_cache = list(original_read_counts.keys())
    
    cache_key = get_cache_key(file_names_for_cache, rarefaction_depth, use_rarefaction)
    
    # Initialize session state cache if needed
    if 'processed_files_cache' not in st.session_state:
        st.session_state.processed_files_cache = {}
    
    # Check cache first
    if cache_key in st.session_state.processed_files_cache:
        file_data = st.session_state.processed_files_cache[cache_key]
    else:
        # Only show processing message when actually processing
        with st.spinner("🔄 Processing files..."):
            file_data = process_files_fresh(uploaded_files, use_rarefaction, rarefaction_depth, original_read_counts)
            st.session_state.processed_files_cache[cache_key] = file_data
            
            # Limit cache size
            if len(st.session_state.processed_files_cache) > 10:
                oldest_key = next(iter(st.session_state.processed_files_cache))
                del st.session_state.processed_files_cache[oldest_key]
    
    if not file_data:
        st.error("No valid insertion files could be loaded.")
        return
    
    # Get selected gene information
    gene_info = annotations[annotations['gene'] == selected_gene]
    
    if gene_info.empty:
        st.error(f"Gene {selected_gene} not found in annotations")
        return
    
    # Display gene information
    st.subheader(f"Gene Information: {selected_gene}")
    gene_data = gene_info.iloc[0]
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Chromosome", gene_data['Chromosome'])
    with col2:
        st.metric("Start", f"{gene_data['Start']:,}")
    with col3:
        st.metric("End", f"{gene_data['End']:,}")
    with col4:
        st.metric("Strand", gene_data['Strand'])
    
    # Use optimized batch processing for overlapping insertions and visualization
    st.subheader("Gene Visualization (All Files)")
    
    # Find overlapping insertions for all files at once - much faster
    overlapping_data = find_overlapping_insertions_batch(gene_info, file_data, chr_mapping)
    
    # Create combined visualization - much faster than merging individual figures
    combined_fig, file_stats_basic = create_combined_gene_visualization(gene_info, overlapping_data)
    
    # Enhance file statistics with additional data
    file_stats = []
    for i, (file_name, insertions_df) in enumerate(sorted(file_data.items())):
        display_name = clean_filename(file_name)
        
        # Get original file name for lookup
        original_file_name = file_name
        if '_rarefied_' in file_name:
            original_file_name = file_name.split('_rarefied_')[0]
            if not uploaded_files:
                original_file_name += '.gz'
        
        stats = {
            'File': display_name,
            'Total Insertions': len(insertions_df),
            'Total Reads': int(insertions_df['Abundance'].sum()),
            'Overlapping Insertions': len(overlapping_data.get(file_name, pd.DataFrame()))
        }
        
        # Add original read count if different from current
        if original_file_name in original_read_counts:
            original_reads = original_read_counts[original_file_name]
            if use_rarefaction and original_reads != stats['Total Reads']:
                stats['Original Total Reads'] = original_reads
        
        file_stats.append(stats)
    
    # Display file statistics
    st.subheader("File Statistics")
    stats_df = pd.DataFrame(file_stats)
    st.dataframe(stats_df, use_container_width=True)
    
    # Display the visualization
    st.plotly_chart(combined_fig, use_container_width=True)
    
    # Show detailed data for each file
    for file_name, overlapping_insertions in overlapping_data.items():
        if not overlapping_insertions.empty:
            display_name = clean_filename(file_name)
            with st.expander(f"Detailed insertions for {display_name} ({len(overlapping_insertions)} insertions)"):
                st.dataframe(overlapping_insertions.sort_values('Abundance', ascending=False))

app()