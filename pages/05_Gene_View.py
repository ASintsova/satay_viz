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

def find_overlapping_insertions(gene_info, insertions_df, chr_mapping, padding=500):
    """Find insertions that overlap with the gene of interest (with padding)"""
    if gene_info.empty or insertions_df.empty:
        return pd.DataFrame()
    
    gene_chr = gene_info['Chromosome'].iloc[0]
    gene_start = gene_info['Start'].iloc[0] - padding
    gene_end = gene_info['End'].iloc[0] + padding
    
    # Map chromosome name to RefSeq ID
    refseq_chr = chr_mapping.get(gene_chr)
    if refseq_chr is None:
        st.warning(f"No RefSeq mapping found for chromosome {gene_chr}")
        return pd.DataFrame()
    
    # Filter insertions on same chromosome (using RefSeq ID)
    chr_insertions = insertions_df[insertions_df['Chromosome'] == refseq_chr]
    
    if chr_insertions.empty:
        return pd.DataFrame()
    
    # Find overlapping insertions (with padding)
    overlapping = chr_insertions[
        (chr_insertions['Start'] <= gene_end) & 
        (chr_insertions['End'] >= gene_start)
    ].copy()
    
    return overlapping

def create_gene_visualization(gene_info, overlapping_insertions, file_name="", y_offset=0, padding=500):
    """Create visualization showing gene rectangle with insertion lines"""
    if gene_info.empty:
        return go.Figure(), 0, 0
    
    gene_start = gene_info['Start'].iloc[0]
    gene_end = gene_info['End'].iloc[0]
    gene_name = gene_info['gene'].iloc[0]
    gene_strand = gene_info['Strand'].iloc[0]
    
    # Expanded region coordinates
    region_start = gene_start - padding
    region_end = gene_end + padding
    
    fig = go.Figure()
    
    # Draw expanded region background
    fig.add_shape(
        type="rect",
        x0=region_start, y0=y_offset-0.8, x1=region_end, y1=y_offset+0.8,
        line=dict(color="lightblue", width=1),
        fillcolor="lightblue",
        opacity=0.2
    )
    
    # Draw gene rectangle
    fig.add_shape(
        type="rect",
        x0=gene_start, y0=y_offset-0.5, x1=gene_end, y1=y_offset+0.5,
        line=dict(color="black", width=2),
        fillcolor="white",
        opacity=1.0
    )
    
    
    # Add file name annotation
    if file_name:
        fig.add_annotation(
            x=region_start,
            y=y_offset+0.6,
            text=file_name,
            showarrow=False,
            font=dict(size=8, color="darkblue"),
            xanchor="left"
        )
    
    # Add strand direction indicator
    if gene_strand == '+':
        # Positive strand: 5' on the left
        fig.add_annotation(
            x=region_start,
            y=y_offset-0.6,
            text="5'",
            showarrow=False,
            font=dict(size=10, color="black"),
            xanchor="left"
        )
    elif gene_strand == '-':
        # Negative strand: 5' on the right
        fig.add_annotation(
            x=region_end,
            y=y_offset-0.6,
            text="5'",
            showarrow=False,
            font=dict(size=10, color="black"),
            xanchor="right"
        )
    
    min_log_abundance = 0
    max_log_abundance = 0
    
    if not overlapping_insertions.empty:
        # Log-scale abundance for color mapping
        log_abundance = np.log10(overlapping_insertions['Abundance'] + 1)  # +1 to handle zeros
        min_log_abundance = log_abundance.min()
        max_log_abundance = log_abundance.max()
        
        # Add insertion lines
        for _, insertion in overlapping_insertions.iterrows():
            # Position lines within the gene rectangle
            x_pos = (insertion['Start'] + insertion['End']) / 2
            
            # Color based on log abundance (red = high abundance, blue = low abundance)
            log_val = np.log10(insertion['Abundance'] + 1)
            if max_log_abundance > min_log_abundance:
                normalized_abundance = (log_val - min_log_abundance) / (max_log_abundance - min_log_abundance)
            else:
                normalized_abundance = 0.5
            
            # Map to red-blue: 0 (blue) to 1 (red)
            red_value = int(normalized_abundance * 255)
            blue_value = int((1 - normalized_abundance) * 255)
            color = f"rgb({red_value},0,{blue_value})"
            
            # Draw vertical line
            fig.add_shape(
                type="line",
                x0=x_pos, y0=y_offset-0.4, x1=x_pos, y1=y_offset+0.4,
                line=dict(color=color, width=2),
                opacity=1.0
            )
            
            # Add invisible scatter point for hover info
            fig.add_trace(go.Scatter(
                x=[x_pos],
                y=[y_offset],
                mode='markers',
                marker=dict(size=0.1, opacity=0),
                hovertemplate=f"File: {file_name}<br>" +
                             f"Position: {insertion['Start']}-{insertion['End']}<br>" +
                             f"Abundance: {insertion['Abundance']}<br>" +
                             f"Log Abundance: {log_val:.2f}<br>" +
                             f"Strand: {insertion['Strand']}<extra></extra>",
                showlegend=False
            ))
    
    # Set x-axis range to show expanded region
    fig.update_xaxes(range=[region_start, region_end])
    
    return fig, min_log_abundance, max_log_abundance

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

def apply_rarefaction_to_dataframe(df, depth, seed=42):
    """Apply rarefaction to a dataframe and return processed dataframe"""
    from rarefy import rarefy
    
    if df.empty:
        return df
    
    try:
        # Apply rarefaction
        rarefied_abundance = rarefy(df['Abundance'].values, depth=depth, seed=seed)
        df = df.copy()
        df['Abundance'] = rarefied_abundance
        
        # Remove zero abundance entries and round to integers
        df = df[df['Abundance'] > 0]
        df['Abundance'] = df['Abundance'].round().astype(int)
        
        return df
    except Exception as e:
        st.error(f"Error during rarefaction: {e}")
        return df

def load_and_process_files(uploaded_files, use_rarefaction, rarefaction_depth):
    """Load files and apply rarefaction if requested"""
    file_data = {}
    original_read_counts = {}
    
    if not uploaded_files:
        # Use rarefied example files
        original_files = [
            "examples/20190221.A-2_noaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz",
            "examples/20190221.A-2_4nMaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz"
        ]
        
        for file_path in original_files:
            file_name = file_path.split('/')[-1]
            insertions_df = load_insertion_file(file_path)
            
            if insertions_df is not None:
                # Store original read count
                original_read_counts[file_name] = int(insertions_df['Abundance'].sum())
                
                if use_rarefaction:
                    insertions_df = apply_rarefaction_to_dataframe(insertions_df, rarefaction_depth)
                    processed_name = f"{file_name}_rarefied_{rarefaction_depth}"
                else:
                    processed_name = file_name
                
                file_data[processed_name] = insertions_df
        
        # Update info message
        if use_rarefaction:
            st.info(f"No files uploaded. Using example files rarefied to {rarefaction_depth:,} reads.")
        else:
            st.info("No files uploaded. Using original example files (no rarefaction).")
    
    else:
        # Load uploaded files
        for uploaded_file in uploaded_files:
            # Get original read count
            original_read_counts[uploaded_file.name] = get_original_read_count(uploaded_file, is_uploaded=True)
            
            # Load and process file
            with tempfile.NamedTemporaryFile(delete=False, suffix=f".{uploaded_file.name.split('.')[-1]}") as tmp_file:
                tmp_file.write(uploaded_file.getvalue())
                tmp_file_path = tmp_file.name
            
            insertions_df = load_insertion_file(tmp_file_path)
            os.unlink(tmp_file_path)
            
            if insertions_df is not None:
                if use_rarefaction:
                    insertions_df = apply_rarefaction_to_dataframe(insertions_df, rarefaction_depth)
                    processed_name = f"{uploaded_file.name}_rarefied_{rarefaction_depth}"
                else:
                    processed_name = uploaded_file.name
                
                file_data[processed_name] = insertions_df
    
    return file_data, original_read_counts

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
    
    # Calculate min reads across files to set rarefaction limit
    min_total_reads = 1000000  # Default high value
    file_read_counts = []
    
    if uploaded_files:
        st.sidebar.subheader("Original File Information")
        for uploaded_file in uploaded_files:
            total_reads = get_original_read_count(uploaded_file, is_uploaded=True)
            file_read_counts.append(total_reads)
            st.sidebar.text(f"{uploaded_file.name}: {total_reads:,} reads")
        min_total_reads = min(file_read_counts) if file_read_counts else 1000000
    else:
        # Use rarefied example files to determine min reads
        example_files = [
            "examples/20190221.A-2_noaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz",
            "examples/20190221.A-2_4nMaF.bam.bed.insertions.sorted.merged.filtered_rarefied_100k.tsv.gz"
        ]
        for file_path in example_files:
            total_reads = get_original_read_count(file_path, is_uploaded=False)
            file_read_counts.append(total_reads)
        min_total_reads = min(file_read_counts) if file_read_counts else 1000000
    
    st.sidebar.header("Rarefaction Settings")
    use_rarefaction = st.sidebar.checkbox("Apply rarefaction", value=True)
    
    # Set max rarefaction depth to minimum total reads across files
    max_rarefaction = max(1000, min_total_reads)  # Ensure at least 1000
    default_rarefaction = min(100000, max_rarefaction)  # Default to 100k or max available
    
    rarefaction_depth = st.sidebar.number_input(
        f"Rarefaction depth (max: {max_rarefaction:,})", 
        min_value=1000, 
        max_value=max_rarefaction, 
        value=default_rarefaction, 
        step=10000
    )
    
    if min_total_reads < 1000000:
        st.sidebar.info(f"Max rarefaction limited to {min_total_reads:,} reads (smallest file)")
    
    # Load and process files
    file_data, original_read_counts = load_and_process_files(uploaded_files, use_rarefaction, rarefaction_depth)
    
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
    
    # Process each file and collect data for combined visualization
    all_figures = []
    all_min_log = []
    all_max_log = []
    file_stats = []
    
    # Sort files by filename for consistent display order
    sorted_file_data = dict(sorted(file_data.items()))
    
    for i, (file_name, insertions_df) in enumerate(sorted_file_data.items()):
        # Find overlapping insertions
        overlapping_insertions = find_overlapping_insertions(gene_info, insertions_df, chr_mapping)
        
        # Clean up file name for display
        display_name = clean_filename(file_name)
        
        # Create visualization for this file
        fig, min_log, max_log = create_gene_visualization(gene_info, overlapping_insertions, display_name, y_offset=i*2)
        all_figures.append(fig)
        all_min_log.append(min_log)
        all_max_log.append(max_log)
        
        # Get original file name for lookup
        original_file_name = file_name
        if '_rarefied_' in file_name:
            # Extract original filename
            original_file_name = file_name.split('_rarefied_')[0]
            if not uploaded_files:
                original_file_name += '.gz'  # Add .gz for example files
        
        # Collect stats
        stats = {
            'File': display_name,
            'Total Insertions': len(insertions_df),
            'Total Reads': int(insertions_df['Abundance'].sum()),
            'Overlapping Insertions': len(overlapping_insertions)
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
    
    # Combine all figures into one
    if all_figures:
        st.subheader("Gene Visualization (All Files)")
        combined_fig = go.Figure()
        
        # Add all traces and shapes from individual figures
        for fig in all_figures:
            for trace in fig.data:
                combined_fig.add_trace(trace)
            for shape in fig.layout.shapes:
                combined_fig.add_shape(shape)
            for annotation in fig.layout.annotations:
                combined_fig.add_annotation(annotation)
        
        # Set up combined layout
        gene_start = gene_info['Start'].iloc[0]
        gene_end = gene_info['End'].iloc[0]
        region_start = gene_start - 500
        region_end = gene_end + 500
        num_files = len(file_data)
        
        # Calculate y-axis range based on actual y_offset values
        # Files are positioned at y_offset = i*2, with gene rectangles spanning ±0.8
        # So range should be from -0.8-padding to (num_files-1)*2+0.8+padding
        y_min = -0.8 - 0.5  # bottom of first file with padding
        y_max = (num_files-1)*2 + 0.8 + 0.5  # top of last file with padding
        
        combined_fig.update_layout(
            title=f"Gene View: {selected_gene} (±500bp)",
            xaxis_title="Genomic Position",
            yaxis_title="Files",
            yaxis=dict(range=[y_min, y_max], showticklabels=False),
            height=200 + num_files * 150,
            showlegend=False,
            template="plotly_white",
            xaxis=dict(range=[region_start, region_end])
        )
        
        # Add colorbar if we have data
        if any(all_min_log) and any(all_max_log):
            global_min_log = min(val for val in all_min_log if val > 0)
            global_max_log = max(all_max_log)
            # Create custom red-blue colorscale (blue = low, red = high)
            red_blue = [[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]
            combined_fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(
                    colorscale=red_blue,
                    cmin=global_min_log,
                    cmax=global_max_log,
                    colorbar=dict(title="Log10(Abundance + 1)")
                ),
                showlegend=False
            ))
        
        st.plotly_chart(combined_fig, use_container_width=True)
        
        # Show detailed data for each file
        for file_name, insertions_df in sorted_file_data.items():
            overlapping_insertions = find_overlapping_insertions(gene_info, insertions_df, chr_mapping)
            if not overlapping_insertions.empty:
                display_name = clean_filename(file_name)
                with st.expander(f"Detailed insertions for {display_name} ({len(overlapping_insertions)} insertions)"):
                    st.dataframe(overlapping_insertions.sort_values('Abundance', ascending=False))

app()