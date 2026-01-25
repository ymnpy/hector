import os, platform, plotly, pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.mesh_utils import free_edges
import warnings

warnings.filterwarnings('ignore')

# ============================================================
# HELPERS
# ============================================================

def my_read_bdf(path):
    """Read BDF file with error handling"""
    try:
        bdf = BDF(debug=False)
        bdf.read_bdf(path, punch=False, xref=False)
    except:
        try: 
            bdf.read_bdf(path, punch=True, xref=False)
        except: 
            bdf.read_bdf(path, punch=False, xref=False)
    return bdf

def plot_free_edges(bdf):
    """Plot lines connecting free edges in the mesh."""
    free_nodes = free_edges.non_paired_edges(bdf)
    xs, ys, zs = [], [], []
    
    for node_pair in free_nodes:
        for nid in node_pair:
            node_coords = bdf.nodes[nid].xyz
            xs.append(node_coords[0])
            ys.append(node_coords[1])
            zs.append(node_coords[2])
        xs.append(None)
        ys.append(None)
        zs.append(None)
    
    return go.Scatter3d(x=xs, y=ys, z=zs, mode='lines',
                       line=dict(color='rgba(0, 0, 0, 0.6)', width=2),
                       showlegend=False, hoverinfo='skip',
                       name='free_edges')


def extract_mesh_data(df, bdf):
    """Extract mesh data from BDF and Excel."""
    
    coords = {'PID': [], 'ELID': [], 'x': [], 'y': [], 'z': [],
              'cx': [], 'cy': [], 'cz': []}
    
    global pid_centroids
    pid_centroids = {}
    
    edges_xyz = {'x': [], 'y': [], 'z': []}
    faces = {'i': [], 'j': [], 'k': []}
    hover_texts = []
    offset = 0
    
    # Define excluded columns
    excluded_columns = ['Property ID', 'ICID', 'Subject_Failure(R)', 'Optimized_Failure(R)']
    intensity_columns = [col for col in df.columns if col not in excluded_columns]
    intensity_values = {col: [] for col in intensity_columns}
    
    # Track which values are non-numeric (for gray coloring)
    non_numeric_mask = {col: [] for col in intensity_columns}
    
    # Create property lookup dictionary (convert once)
    property_data = {}
    for idx, row in df.iterrows():
        try:
            pid = int(float(row['Property ID']))
            property_data[pid] = row
        except (ValueError, KeyError) as e:
            print(f"Warning: Skipping row {idx} - invalid Property ID: {e}")
            continue
    
    elements_processed = 0
    
    for eid, element in bdf.elements.items():
        pid = element.pid
        
        # Skip if property not in Excel data
        if pid not in property_data:
            continue
        
        elements_processed += 1
        
        coords['PID'].append(pid)
        coords['ELID'].append(eid)
        
        # Get property data
        rf_data = property_data[pid]
        
        # Build hover info
        hover_info = [f"<b>Element ID:</b> {eid}", f"<b>Property ID:</b> {pid}"]
        for col in intensity_columns:
            try:
                value = float(rf_data[col])
                hover_info.append(f'<b>{col}:</b> {value:.4f}')
            except (ValueError, TypeError):
                hover_info.append(f'<b>{col}:</b> {rf_data[col]}')
        
        hover_text = "<br>".join(hover_info)
        
        # Get element type and nodes
        element_type = element.type
        element_coords = []
        
        # Collect node coordinates
        for nid in element.nodes:
            node_xyz = bdf.nodes[nid].xyz
            coords['x'].append(node_xyz[0])
            coords['y'].append(node_xyz[1])
            coords['z'].append(node_xyz[2])
            element_coords.append(node_xyz)
        
        # Calculate centroid
        centroid = np.mean(element_coords, axis=0)
        coords['cx'].append(centroid[0])
        coords['cy'].append(centroid[1])
        coords['cz'].append(centroid[2])
        pid_centroids.setdefault(pid, []).append(centroid)
        
        # Process based on element type
        if element_type == 'CTRIA3':
            # Store intensity values - ONE per element
            for col in intensity_columns:
                try:
                    value = float(rf_data[col])
                    intensity_values[col].append(value)
                    non_numeric_mask[col].append(False)
                except (ValueError, TypeError):
                    # Use a placeholder value (will be colored gray)
                    intensity_values[col].append(-999999)  # Sentinel value
                    non_numeric_mask[col].append(True)
            
            hover_texts.append(hover_text)
            
            # Add triangle edges
            for edge in [[offset, offset+1], [offset+1, offset+2], [offset+2, offset]]:
                edges_xyz['x'].extend([coords['x'][edge[0]], coords['x'][edge[1]], None])
                edges_xyz['y'].extend([coords['y'][edge[0]], coords['y'][edge[1]], None])
                edges_xyz['z'].extend([coords['z'][edge[0]], coords['z'][edge[1]], None])
            
            # Add triangle face
            faces['i'].append(offset)
            faces['j'].append(offset+1)
            faces['k'].append(offset+2)
            offset += 3
            
        elif element_type == 'CQUAD4':
            # CQUAD4 creates TWO faces, so duplicate intensity values
            for col in intensity_columns:
                try:
                    val = float(rf_data[col])
                    intensity_values[col].extend([val, val])
                    non_numeric_mask[col].extend([False, False])
                except (ValueError, TypeError):
                    intensity_values[col].extend([-999999, -999999])
                    non_numeric_mask[col].extend([True, True])
            
            # Duplicate hover text for both faces
            hover_texts.extend([hover_text, hover_text])
            
            # Add quad edges
            for edge in [[offset, offset+1], [offset+1, offset+2], 
                        [offset+2, offset+3], [offset+3, offset]]:
                edges_xyz['x'].extend([coords['x'][edge[0]], coords['x'][edge[1]], None])
                edges_xyz['y'].extend([coords['y'][edge[0]], coords['y'][edge[1]], None])
                edges_xyz['z'].extend([coords['z'][edge[0]], coords['z'][edge[1]], None])
            
            # Add quad faces (two triangles)
            faces['i'].extend([offset, offset+2])
            faces['j'].extend([offset+1, offset+3])
            faces['k'].extend([offset+2, offset])
            offset += 4
    
    print(f"Processed {elements_processed} elements, created {len(faces['i'])} faces")
    
    return coords, faces, edges_xyz, intensity_values, hover_texts, intensity_columns, non_numeric_mask


def create_property_labels(coords, intensity_values, intensity_columns):
    """Create label traces for each property at centroid locations - showing intensity values."""
    
    label_traces = {}
    
    # Get unique PIDs and their centroids
    unique_pids = list(set(coords['PID']))
    
    for col in intensity_columns:
        label_x, label_y, label_z, label_text = [], [], [], []
        
        for pid in unique_pids:
            # Find all indices for this PID
            pid_indices = [i for i, p in enumerate(coords['PID']) if p == pid]
            
            if not pid_indices:
                continue
            
            # Get centroid (average of all centroids for this PID)
            cx = np.mean([coords['cx'][i] for i in pid_indices])
            cy = np.mean([coords['cy'][i] for i in pid_indices])
            cz = np.mean([coords['cz'][i] for i in pid_indices])
            
            # Get value for this property from THIS COLUMN
            if pid_indices[0] < len(intensity_values[col]):
                value = intensity_values[col][pid_indices[0]]
                
                label_x.append(cx)
                label_y.append(cy)
                label_z.append(cz)
                
                # Show the actual intensity value, or "N/A" for non-numeric
                if value == -999999:
                    label_text.append("N/A")
                else:
                    label_text.append(f"{value:.2f}")
        
        # Create scatter trace for labels
        label_trace = go.Scatter3d(
            x=label_x, y=label_y, z=label_z,
            mode='text',
            text=label_text,
            textfont=dict(
                size=11, 
                color='white',
                family='Arial Black, sans-serif'
            ),
            textposition='middle center',
            showlegend=False,
            hoverinfo='skip',
            name=f'labels_{col}',
            visible=False  # Hidden by default
        )
        
        label_traces[col] = label_trace
    
    return label_traces


def get_discrete_colorscale(values):
    """Generate discrete colorscale for categorical/ID data"""
    unique_vals = sorted(list(set(values)))
    n_colors = len(unique_vals)
    
    # Use vibrant, distinguishable colors
    if n_colors <= 10:
        colors = px.colors.qualitative.Bold
    elif n_colors <= 24:
        colors = px.colors.qualitative.Dark24
    else:
        # For many values, use a continuous colorscale but make it vibrant
        colors = px.colors.sequential.Rainbow
    
    # Create discrete colorscale
    colorscale = []
    for i, val in enumerate(unique_vals):
        color = colors[i % len(colors)]
        position = i / max(n_colors - 1, 1)
        colorscale.append((position, color))
    
    return colorscale, unique_vals


def create_mesh_visualization(df, coords, faces, edges_xyz, intensity_values, hover_texts, intensity_columns, non_numeric_mask):
    """Create mesh plot with dynamically generated intensity meshes and buttons."""
    
    mesh_traces = []
    
    # Remove unwanted columns from visualization
    removed_ones = ["Subcase3D", "Subcase_2D", "LCID"]
    intensity_columns = [col for col in intensity_columns if col not in removed_ones]
    
    # Define color scale categories
    discrete_ones = ['C_(kg)_default', 'C_(kg)_optimized', 'mass_(kg)_default', 'mass_(kg)_optimizer','min_RF']
    reversed_ones = ['Property_ID', 'ICID', 'UCID', 'V_min_size', 'V_mid_size', 'V_max_size']
    
    # Define which columns should use discrete coloring
    discrete_id_columns = ['Subcase_ID', 'UCID', 'Property_ID']
    discrete_rf_columns = ['min_RF', 'RF_combined', 'RF_bending', 'RF_membrane']
    
    # Create a mesh for each intensity column
    for i, col in enumerate(intensity_columns):
        intensity_data = intensity_values.get(col, [])
        mask = non_numeric_mask.get(col, [])
        
        # Separate numeric and non-numeric data
        numeric_data = [v for v, is_non_numeric in zip(intensity_data, mask) if not is_non_numeric]
        
        # Check if this column should use discrete coloring
        use_discrete = any(keyword in col for keyword in discrete_id_columns + discrete_rf_columns)
        
        if use_discrete and numeric_data:
            # Use discrete colorscale for IDs and RF values
            if any(keyword in col for keyword in discrete_id_columns):
                colorscale, unique_vals = get_discrete_colorscale(numeric_data)
                colorbar_config = dict(
                    title=dict(text=col, side='right', font=dict(size=12, color='#333')),
                    thickness=20,
                    len=0.85,
                    x=1.02,
                    xanchor='left',
                    tickmode='array',
                    tickvals=unique_vals,
                    ticktext=[str(int(v)) for v in unique_vals],
                    tickfont=dict(size=9)
                )
            else:
                # For RF values - STEPPED discrete rainbow colorscale with ranges
                colorscale = [
                    (0.0, 'grey'),      # Range: 0.0 - 0.125
                    (0.125, 'grey'),
                    (0.125, 'blue'),    # Range: 0.125 - 0.25
                    (0.25, 'blue'),
                    (0.25, 'cyan'),     # Range: 0.25 - 0.375
                    (0.375, 'cyan'),
                    (0.375, 'green'),   # Range: 0.375 - 0.5
                    (0.5, 'green'),
                    (0.5, 'yellow'),    # Range: 0.5 - 0.625
                    (0.625, 'yellow'),
                    (0.625, 'orange'),  # Range: 0.625 - 0.75
                    (0.75, 'orange'),
                    (0.75, 'magenta'),  # Range: 0.75 - 0.875
                    (0.875, 'magenta'),
                    (0.875, 'red'),     # Range: 0.875 - 1.0
                    (1.0, 'red')
                ]
                data_min, data_max = min(numeric_data), max(numeric_data)
                # Create 8 discrete ranges
                n_ranges = 8
                tick_vals = [data_min + (data_max - data_min) * j / n_ranges for j in range(n_ranges + 1)]
                colorbar_config = dict(
                    title=dict(text=col, side='right', font=dict(size=12, color='#333')),
                    thickness=20,
                    len=0.85,
                    x=1.02,
                    xanchor='left',
                    tickmode='array',
                    tickvals=tick_vals,
                    ticktext=[f'{v:.2f}' for v in tick_vals],
                    tickfont=dict(size=9)
                )
        elif col in reversed_ones:
            # STEPPED discrete rainbow colorscale with ranges
            colorscale = [
                (0.0, 'grey'),      # Range: 0.0 - 0.125
                (0.125, 'grey'),
                (0.125, 'blue'),    # Range: 0.125 - 0.25
                (0.25, 'blue'),
                (0.25, 'cyan'),     # Range: 0.25 - 0.375
                (0.375, 'cyan'),
                (0.375, 'green'),   # Range: 0.375 - 0.5
                (0.5, 'green'),
                (0.5, 'yellow'),    # Range: 0.5 - 0.625
                (0.625, 'yellow'),
                (0.625, 'orange'),  # Range: 0.625 - 0.75
                (0.75, 'orange'),
                (0.75, 'magenta'),  # Range: 0.75 - 0.875
                (0.875, 'magenta'),
                (0.875, 'red'),     # Range: 0.875 - 1.0
                (1.0, 'red')
            ]
            if numeric_data:
                data_min, data_max = min(numeric_data), max(numeric_data)
                n_ranges = 8
                tick_vals = [data_min + (data_max - data_min) * j / n_ranges for j in range(n_ranges + 1)]
                colorbar_config = dict(
                    title=dict(text=col, side='right', font=dict(size=12, color='#333')),
                    thickness=20,
                    len=0.85,
                    x=1.02,
                    xanchor='left',
                    tickmode='array',
                    tickvals=tick_vals,
                    ticktext=[f'{v:.2f}' for v in tick_vals],
                    tickfont=dict(size=9)
                )
            else:
                colorbar_config = dict(title=dict(text=col, side='right'))
        elif col in discrete_ones:
            colorscale = [
                (0, "#37474f"),
                (0.15, "#1e88e5"),
                (0.35, "#00acc1"),
                (0.5, "#7cb342"),
                (0.65, "#fdd835"),
                (0.8, "#fb8c00"),
                (1.0, "#e53935")
            ]
            if numeric_data:
                data_min, data_max = min(numeric_data), max(numeric_data)
                colorbar_config = dict(
                    title=dict(text=col, side='right', font=dict(size=12, color='#333')),
                    thickness=20,
                    len=0.85,
                    x=1.02,
                    xanchor='left',
                    tickvals=[data_min, (data_min + data_max) / 2, data_max],
                    ticktext=[f'{data_min:.2f}', f'{(data_min + data_max) / 2:.2f}', f'{data_max:.2f}'],
                    tickfont=dict(size=10)
                )
            else:
                colorbar_config = dict(title=dict(text=col, side='right'))
        else:
            # Modern stress/default colorscale
            colorscale = [
                (0, "#263238"),
                (0.5, "#43a047"),
                (0.75, "#fdd835"),
                (0.9, "#fb8c00"),
                (1, "#d32f2f")
            ]
            if numeric_data:
                data_min, data_max = min(numeric_data), max(numeric_data)
                colorbar_config = dict(
                    title=dict(text=col, side='right', font=dict(size=12, color='#333')),
                    thickness=20,
                    len=0.85,
                    x=1.02,
                    xanchor='left',
                    tickvals=[data_min, (data_min + data_max) / 2, data_max],
                    ticktext=[f'{data_min:.2f}', f'{(data_min + data_max) / 2:.2f}', f'{data_max:.2f}'],
                    tickfont=dict(size=10)
                )
            else:
                colorbar_config = dict(title=dict(text=col, side='right'))
        
        # Replace non-numeric values with min value for gray coloring
        if numeric_data:
            min_val = min(numeric_data)
            processed_data = [min_val if is_non_numeric else v 
                            for v, is_non_numeric in zip(intensity_data, mask)]
        else:
            processed_data = intensity_data
        
        if numeric_data:
            print(f"Column '{col}': {len(numeric_data)} numeric values, {sum(mask)} non-numeric (gray), range: {min(numeric_data):.4f} to {max(numeric_data):.4f}")
        
        mesh = go.Mesh3d(
            x=coords['x'], 
            y=coords['y'], 
            z=coords['z'],
            i=faces['i'], 
            j=faces['j'], 
            k=faces['k'],
            intensity=processed_data,
            colorscale=colorscale,
            intensitymode='cell',
            showscale=True,
            colorbar=colorbar_config,
            hovertemplate='%{text}<extra></extra>',
            text=hover_texts,
            name=f"mesh_{col}",
            visible=(i == 0),
            lighting=dict(ambient=0.7, diffuse=0.9, specular=0.5, roughness=0.3, fresnel=0.2),
            flatshading=False,
            opacity=1.0
        )
        mesh_traces.append(mesh)
    
    # Create edges trace with modern styling
    edges_trace = go.Scatter3d(
        x=edges_xyz['x'], y=edges_xyz['y'], z=edges_xyz['z'],
        mode='lines',
        line=dict(color='rgba(50, 50, 50, 0.3)', width=0.8),
        showlegend=False,
        hoverinfo='skip',
        name='edge_traces'
    )
    
    # Create property labels
    label_traces = create_property_labels(coords, intensity_values, intensity_columns)
    
    # Create buttons dynamically for each intensity column
    buttons = []
    num_mesh_traces = len(intensity_columns)
    
    for i, col in enumerate(intensity_columns):
        # Visibility: [meshes..., edges, free_edges, labels...]
        visible = [False] * num_mesh_traces  # All meshes
        visible[i] = True  # Current mesh
        visible.extend([True, True])  # edges_trace and free_edges_trace
        # Add label visibility - only show label for current column
        label_visible = [j == i for j in range(len(intensity_columns))]
        visible.extend(label_visible)
        
        buttons.append({
            'label': col,
            'method': 'update',
            'args': [{'visible': visible}, {'title': None}]
        })
    
    updatemenus = [{
        'buttons': buttons,
        'direction': 'down',
        'showactive': True,
        'active': 0,
        'x': 0.87,
        'xanchor': "left",
        'y': 0.98,
        'yanchor': "top",
        'bgcolor': 'rgba(255, 255, 255, 0.9)',
        'bordercolor': '#1976d2',
        'borderwidth': 2,
        'font': dict(size=11, color='#333', family='Arial, sans-serif')
    }]
    
    return mesh_traces, edges_trace, label_traces, updatemenus, intensity_columns


def main(bdf_file, excel_file):
    """Load YTE summary excel and plots 3D mesh plot w/ plotly by also using pyNastran"""
    
    # Load data
    df = pd.read_excel(excel_file)
    bdf_model = my_read_bdf(bdf_file)
    
    # Extract mesh data
    coords, faces, edges_xyz, intensity_values, hover_texts, intensity_columns, non_numeric_mask = extract_mesh_data(
        df, bdf_model
    )
    
    if len(faces['i']) == 0:
        print("ERROR: No faces were created!")
        return None
    
    # Create visualization with dynamic buttons
    mesh_traces, edges_trace, label_traces, updatemenus_main, viz_columns = create_mesh_visualization(
        df, coords, faces, edges_xyz, intensity_values, hover_texts, intensity_columns, non_numeric_mask
    )
    
    # Add free edges
    free_edges_trace = plot_free_edges(bdf_model)
    
    # Create figure
    fig = go.Figure()
    
    # Add all traces in order: meshes, edges, free_edges, labels
    for trace in mesh_traces:
        fig.add_trace(trace)
    fig.add_trace(edges_trace)
    fig.add_trace(free_edges_trace)
    for col in viz_columns:
        fig.add_trace(label_traces[col])
    
    # Get the number of traces
    num_mesh_traces = len(mesh_traces)
    num_total_traces = len(fig.data)
    
    # Additional control buttons - modern styled and aligned
    updatemenus_reset = dict(
        type='buttons',
        active=0,
        showactive=True,
        buttons=[dict(
            args=[{
                'scene.camera.center': {'x': 0, 'y': 0, 'z': 0},
                'scene.camera.eye': {'x': 1.25, 'y': 1.25, 'z': 1.25},
                'scene.camera.up': {'x': 0, 'y': 0, 'z': 1}
            }],
            label='↻ Reset',
            method='relayout'
        )],
        direction='down',
        x=0.01,
        y=0.02,
        xanchor='left',
        yanchor='bottom',
        font=dict(size=11, color='#333', family='Arial, sans-serif'),
        bgcolor='rgba(255, 255, 255, 0.95)',
        bordercolor='#1976d2',
        borderwidth=2
    )
    
    updatemenus_opacity = dict(
        type='buttons',
        showactive=True,
        active=0,
        x=0.095,
        xanchor='left',
        y=0.02,
        yanchor='bottom',
        buttons=[dict(
            label='◐ Opacity',
            method='restyle',
            args=[{'opacity': 0.5}, list(range(num_mesh_traces))],
            args2=[{'opacity': 1}, list(range(num_mesh_traces))]
        )],
        font=dict(color='#333', size=11, family='Arial, sans-serif'),
        bgcolor='rgba(255, 255, 255, 0.95)',
        bordercolor='#1976d2',
        borderwidth=2
    )
    
    updatemenus_edges = dict(
        type='buttons',
        showactive=True,
        active=0,
        x=0.215,
        xanchor='left',
        y=0.02,
        yanchor='bottom',
        buttons=[dict(
            label='◻ Edges',
            method='restyle',
            args=[{'visible': False}, [num_mesh_traces]],  # Hide edges trace
            args2=[{'visible': True}, [num_mesh_traces]]    # Show edges trace
        )],
        font=dict(color='#333', size=11, family='Arial, sans-serif'),
        bgcolor='rgba(255, 255, 255, 0.95)',
        bordercolor='#1976d2',
        borderwidth=2
    )
    
    updatemenus_theme = dict(
        type='buttons',
        direction='right',
        active=0,
        showactive=True,
        x=0.32,
        xanchor='left',
        y=0.02,
        yanchor='bottom',
        buttons=[dict(
            label='☀ Theme',
            method='relayout',
            args=[{'template': 'plotly_white', 'scene.bgcolor': '#f5f5f5'}],
            args2=[{'template': 'plotly_dark', 'scene.bgcolor': '#1a1a1a'}]
        )],
        font=dict(color='#333', size=11, family='Arial, sans-serif'),
        bgcolor='rgba(255, 255, 255, 0.95)',
        bordercolor='#1976d2',
        borderwidth=2
    )
    
    # Labels toggle button
    updatemenus_labels = dict(
        type='buttons',
        showactive=True,
        active=0,
        x=0.43,
        xanchor='left',
        y=0.02,
        yanchor='bottom',
        buttons=[dict(
            label='🏷 Labels',
            method='restyle',
            args=[{'visible': True}, list(range(num_mesh_traces + 2, num_total_traces))],
            args2=[{'visible': False}, list(range(num_mesh_traces + 2, num_total_traces))]
        )],
        font=dict(color='#333', size=11, family='Arial, sans-serif'),
        bgcolor='rgba(255, 255, 255, 0.95)',
        bordercolor='#1976d2',
        borderwidth=2
    )
    
    # Configure layout with modern styling
    fig.update_layout(
        scene=dict(
            aspectmode='data',
            xaxis=dict(visible=False, showgrid=False),
            yaxis=dict(visible=False, showgrid=False),
            zaxis=dict(visible=False, showgrid=False),
            bgcolor='#f5f5f5'
        ),
        template='plotly_white',
        updatemenus=[updatemenus_main[0], updatemenus_reset, updatemenus_opacity, updatemenus_edges, updatemenus_theme, updatemenus_labels],
        hovermode='closest',
        title=None,
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0),
        autosize=True,
        paper_bgcolor='#fafafa',
        font=dict(family='Arial, sans-serif', color='#333')
    )
    
    output_file = 'plotly_output.html'
    fig.write_html(output_file)
    print(f"Complete: plot written to {output_file}")
    
    return fig


if __name__ == "__main__":
    path1 = r'C:\Users\User\Desktop\vaeridion\HM\updated_model_220619_2249\updated_statics.bdf'
    path2 = r'C:/Users/User/Desktop/VSCODE/13_PLOTLY/10_minRF_plot/New folder/final.xlsx'
    
    main(path1, path2)
