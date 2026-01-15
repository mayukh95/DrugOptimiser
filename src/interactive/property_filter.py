"""
Interactive Property Filter with histograms and sliders.
"""

import ipywidgets as widgets
from IPython.display import display, clear_output

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


class InteractivePropertyFilter:
    """Filter molecules by properties with live histogram visualization."""
    
    def __init__(self, molecules_list=None):
        """
        Initialize property filter.
        
        Parameters:
        -----------
        molecules_list : list
            List of DrugMolecule objects
        """
        self.molecules = molecules_list if molecules_list else []
        self.filtered_molecules = list(self.molecules)
        self._setup_widgets()
    
    def _setup_widgets(self):
        """Create filter widgets."""
        self.mw_range = widgets.FloatRangeSlider(
            value=[0, 1000],
            min=0, max=1500,
            step=10,
            description='MW:',
            style={'description_width': '80px'},
            layout=widgets.Layout(width='350px'),
            readout_format='.0f'
        )
        
        self.logp_range = widgets.FloatRangeSlider(
            value=[-5, 10],
            min=-10, max=15,
            step=0.5,
            description='LogP:',
            style={'description_width': '80px'},
            layout=widgets.Layout(width='350px')
        )
        
        self.tpsa_range = widgets.FloatRangeSlider(
            value=[0, 200],
            min=0, max=300,
            step=5,
            description='TPSA:',
            style={'description_width': '80px'},
            layout=widgets.Layout(width='350px')
        )
        
        self.qed_range = widgets.FloatRangeSlider(
            value=[0, 1],
            min=0, max=1,
            step=0.05,
            description='QED:',
            style={'description_width': '80px'},
            layout=widgets.Layout(width='350px')
        )
        
        self.hbd_max = widgets.IntSlider(
            value=10,
            min=0, max=20,
            description='Max HBD:',
            style={'description_width': '80px'},
            layout=widgets.Layout(width='350px')
        )
        
        self.hba_max = widgets.IntSlider(
            value=15,
            min=0, max=30,
            description='Max HBA:',
            style={'description_width': '80px'},
            layout=widgets.Layout(width='350px')
        )
        
        self.filter_btn = widgets.Button(
            description='🔍 Apply Filters',
            button_style='primary',
            layout=widgets.Layout(width='150px')
        )
        self.filter_btn.on_click(self._apply_filters)
        
        self.reset_btn = widgets.Button(
            description='🔄 Reset',
            button_style='warning',
            layout=widgets.Layout(width='100px')
        )
        self.reset_btn.on_click(self._reset_filters)
        
        self.result_output = widgets.Output()
        self.histogram_output = widgets.Output()
        self.molecule_list_output = widgets.Output(layout=widgets.Layout(
            max_height='300px', overflow_y='auto'
        ))
    
    def load_molecules(self, molecules_list):
        """
        Load molecules into filter.
        
        Parameters:
        -----------
        molecules_list : list
            List of DrugMolecule objects
        """
        self.molecules = molecules_list
        self.filtered_molecules = list(molecules_list)
        self._apply_filters()
    
    def _apply_filters(self, b=None):
        """Apply current filter settings."""
        self.filtered_molecules = []
        
        for mol in self.molecules:
            if not mol.properties:
                mol.calculate_descriptors()
            
            props = mol.properties
            
            # Check all filter conditions
            if not (self.mw_range.value[0] <= props.get('MW', 0) <= self.mw_range.value[1]):
                continue
            if not (self.logp_range.value[0] <= props.get('LogP', 0) <= self.logp_range.value[1]):
                continue
            if not (self.tpsa_range.value[0] <= props.get('TPSA', 0) <= self.tpsa_range.value[1]):
                continue
            if not (self.qed_range.value[0] <= props.get('QED', 0) <= self.qed_range.value[1]):
                continue
            if props.get('HBD', 0) > self.hbd_max.value:
                continue
            if props.get('HBA', 0) > self.hba_max.value:
                continue
            
            self.filtered_molecules.append(mol)
        
        self._show_results()
    
    def _reset_filters(self, b=None):
        """Reset all filters."""
        self.mw_range.value = [0, 1000]
        self.logp_range.value = [-5, 10]
        self.tpsa_range.value = [0, 200]
        self.qed_range.value = [0, 1]
        self.hbd_max.value = 10
        self.hba_max.value = 15
        self._apply_filters()
    
    def _show_results(self):
        """Display filtered results."""
        with self.result_output:
            clear_output(wait=True)
            
            total = len(self.molecules)
            filtered = len(self.filtered_molecules)
            pct = (filtered / total * 100) if total > 0 else 0
            
            print(f"📊 Filtered Molecules: {filtered}/{total} ({pct:.1f}%)")
            print("─" * 50)
            
            for mol in self.filtered_molecules[:10]:
                props = mol.properties
                print(f"  • {mol.name}: MW={props.get('MW', 0):.1f}, "
                      f"LogP={props.get('LogP', 0):.2f}, QED={props.get('QED', 0):.3f}")
            
            if len(self.filtered_molecules) > 10:
                print(f"  ... and {len(self.filtered_molecules) - 10} more")
        
        # Populate molecule list table
        self._show_molecule_list()
        
        self._plot_histograms()
    
    def _show_molecule_list(self):
        """Display filtered molecules as a styled table."""
        from IPython.display import HTML
        
        with self.molecule_list_output:
            clear_output(wait=True)
            
            if not self.filtered_molecules:
                display(HTML("<p style='color:#888; padding:20px;'>No molecules match the current filters.</p>"))
                return
            
            # Build HTML table
            html = """
            <style>
                .mol-table { width:100%; border-collapse:collapse; font-size:13px; }
                .mol-table th { background:linear-gradient(135deg, #667eea, #764ba2); color:white; 
                               padding:10px; text-align:left; position:sticky; top:0; }
                .mol-table td { padding:8px 10px; border-bottom:1px solid #eee; }
                .mol-table tr:hover { background:#f5f5ff; }
                .mol-table tr:nth-child(even) { background:#fafafa; }
                .pass { color:#10b981; font-weight:600; }
                .fail { color:#ef4444; }
                .qed-high { background:#d1fae5; color:#065f46; padding:2px 6px; border-radius:4px; }
                .qed-med { background:#fef3c7; color:#92400e; padding:2px 6px; border-radius:4px; }
                .qed-low { background:#fee2e2; color:#991b1b; padding:2px 6px; border-radius:4px; }
            </style>
            <table class='mol-table'>
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Molecule</th>
                        <th>MW</th>
                        <th>LogP</th>
                        <th>TPSA</th>
                        <th>HBD</th>
                        <th>HBA</th>
                        <th>QED</th>
                        <th>Lipinski</th>
                    </tr>
                </thead>
                <tbody>
            """
            
            for i, mol in enumerate(self.filtered_molecules, 1):
                props = mol.properties
                mw = props.get('MW', 0)
                logp = props.get('LogP', 0)
                tpsa = props.get('TPSA', 0)
                hbd = props.get('HBD', 0)
                hba = props.get('HBA', 0)
                qed = props.get('QED', 0)
                
                # Lipinski check
                violations = 0
                if mw > 500: violations += 1
                if logp > 5: violations += 1
                if hbd > 5: violations += 1
                if hba > 10: violations += 1
                lipinski = f"<span class='pass'>✓ Pass</span>" if violations == 0 else f"<span class='fail'>✗ {violations}v</span>"
                
                # QED styling
                if qed >= 0.67:
                    qed_class = 'qed-high'
                elif qed >= 0.4:
                    qed_class = 'qed-med'
                else:
                    qed_class = 'qed-low'
                
                html += f"""
                    <tr>
                        <td>{i}</td>
                        <td><b>{mol.name}</b></td>
                        <td>{mw:.1f}</td>
                        <td>{logp:.2f}</td>
                        <td>{tpsa:.1f}</td>
                        <td>{hbd}</td>
                        <td>{hba}</td>
                        <td><span class='{qed_class}'>{qed:.3f}</span></td>
                        <td>{lipinski}</td>
                    </tr>
                """
            
            html += "</tbody></table>"
            display(HTML(html))
    
    def _plot_histograms(self):
        """Plot property distributions."""
        with self.histogram_output:
            clear_output(wait=True)
            
            if not self.filtered_molecules:
                print("No molecules match current filters")
                return
            
            if not PLOTLY_AVAILABLE:
                print("Plotly not available for histograms")
                return
            
            props_to_plot = ['MW', 'LogP', 'QED', 'TPSA']
            values = {p: [] for p in props_to_plot}
            
            for mol in self.filtered_molecules:
                for p in props_to_plot:
                    values[p].append(mol.properties.get(p, 0))
            
            fig = make_subplots(rows=2, cols=2, subplot_titles=props_to_plot)
            
            colors = ['#667eea', '#f56565', '#38a169', '#ed8936']
            
            for i, (prop, color) in enumerate(zip(props_to_plot, colors)):
                row = i // 2 + 1
                col = i % 2 + 1
                
                fig.add_trace(
                    go.Histogram(x=values[prop], nbinsx=15, 
                                marker_color=color, opacity=0.7,
                                name=prop),
                    row=row, col=col
                )
            
            fig.update_layout(
                height=400,
                width=600,
                showlegend=False,
                title_text="Property Distributions",
                margin=dict(l=40, r=40, t=60, b=40)
            )
            
            fig.show()
    
    def get_filtered_molecules(self):
        """Get the currently filtered molecules."""
        return self.filtered_molecules
    
    def display(self):
        """Display the filter interface."""
        header = widgets.HTML("""
        <div style='background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); 
                    padding: 16px; border-radius: 12px; margin-bottom: 16px;'>
            <h3 style='color: white; margin: 0;'>🎛️ Interactive Property Filter</h3>
            <p style='color: rgba(255,255,255,0.85); margin: 8px 0 0 0; font-size: 13px;'>
                Filter molecules by drug-like properties with live visualization
            </p>
        </div>
        """)
        
        filters_left = widgets.VBox([
            self.mw_range,
            self.logp_range,
            self.tpsa_range
        ])
        
        filters_right = widgets.VBox([
            self.qed_range,
            self.hbd_max,
            self.hba_max
        ])
        
        filters_row = widgets.HBox([filters_left, filters_right], 
                                   layout=widgets.Layout(gap='20px'))
        
        buttons_row = widgets.HBox([self.filter_btn, self.reset_btn],
                                   layout=widgets.Layout(gap='12px', margin='16px 0'))
        
        results_row = widgets.HBox([
            widgets.VBox([
                widgets.HTML("<b>📋 Summary</b>"),
                self.result_output
            ], layout=widgets.Layout(width='380px', padding='12px',
                                    border='1px solid #e0e0e0', border_radius='8px')),
            widgets.VBox([
                widgets.HTML("<b>📊 Distributions</b>"),
                self.histogram_output
            ], layout=widgets.Layout(width='650px', padding='12px',
                                    border='1px solid #e0e0e0', border_radius='8px'))
        ], layout=widgets.Layout(gap='16px'))
        
        # Molecule list panel
        molecule_list_panel = widgets.VBox([
            widgets.HTML("""
                <div style='display:flex; align-items:center; gap:8px; margin-bottom:8px;'>
                    <b>📝 Filtered Molecules List</b>
                    <span style='color:#888; font-size:12px;'>(scrollable)</span>
                </div>
            """),
            self.molecule_list_output
        ], layout=widgets.Layout(
            width='100%', 
            padding='12px',
            border='2px solid #667eea', 
            border_radius='8px',
            margin='16px 0 0 0'
        ))
        
        interface = widgets.VBox([
            header,
            filters_row,
            buttons_row,
            results_row,
            molecule_list_panel
        ], layout=widgets.Layout(max_width='1100px'))
        
        # Show initial results
        if self.molecules:
            self._apply_filters()
        
        return interface