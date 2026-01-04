"""
3D Protein Rendering using py3Dmol for web-based visualization
"""

import numpy as np
import tempfile
import os
import io
from Bio.PDB import MMCIFIO
from typing import Dict, Any, List, Tuple, Optional
from ..core.protein_structure import ProteinStructure
from ..core.pigment_system import PigmentSystem


class Protein3DRenderer:
    """
    Renderer for 3D protein visualization using py3Dmol
    """
    
    def __init__(self, structure: ProteinStructure, pigment_system: PigmentSystem):
        """
        Initialize renderer with protein structure and pigment system
        
        Args:
            structure: ProteinStructure object
            pigment_system: PigmentSystem object
        """
        self.structure = structure
        self.pigment_system = pigment_system
        self.temp_files = []
    
    def generate_html_view(self, width: str = "100%", height: str = "100vh") -> str:
        """
        Generate complete HTML view with embedded 3D visualization
        """

        # Prepare structure and pigment data for embedding in the HTML/JS template
        # Use CIF format to handle extended structures and large coordinates
        cif_content = self.generate_cif_content()
        cif_content_escaped = cif_content.replace('\\', '\\\\').replace('\n', '\\n').replace("'", "\\'")
        pigment_info = self.get_pigment_info_js()

        html_template = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Alprotein 3D Viewer</title>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                :root {
                    --bg: #f6f8fb;
                    --panel: #ffffff;
                    --border: #e6e8ef;
                    --text: #1f2a3d;
                    --muted: #6b7280;
                    --accent: #2563eb;
                    --accent-soft: #e8f0ff;
                    --shadow: 0 12px 30px rgba(31, 42, 61, 0.12);
                    --radius: 10px;
                    --font: 'Inter', 'Segoe UI', system-ui, -apple-system, sans-serif;
                }
                * { box-sizing: border-box; }
                body {
                    margin: 0;
                    padding: 0;
                    background: var(--bg);
                    font-family: var(--font);
                    color: var(--text);
                    overflow: hidden;
                }
                #viewer {
                    width: __WIDTH__;
                    height: __HEIGHT__;
                    position: relative;
                }
                #info-panel {
                    position: absolute;
                    top: 16px;
                    right: 16px;
                    background: var(--panel);
                    border: 1px solid var(--border);
                    border-radius: var(--radius);
                    padding: 14px 16px;
                    box-shadow: var(--shadow);
                    z-index: 1000;
                    min-width: 220px;
                }
                #info-panel h3 {
                    margin: 0 0 8px 0;
                    font-size: 15px;
                    font-weight: 700;
                    display: flex;
                    align-items: center;
                    gap: 8px;
                }
                #info-panel .dot {
                    width: 10px;
                    height: 10px;
                    background: var(--accent);
                    border-radius: 50%;
                    box-shadow: 0 0 0 4px var(--accent-soft);
                }
                #info-panel .stat-row {
                    display: flex;
                    gap: 12px;
                    font-size: 12px;
                    color: var(--muted);
                    margin-bottom: 8px;
                }
                #info-panel .stat {
                    display: flex;
                    gap: 4px;
                    align-items: center;
                }
                #info-panel .tip {
                    font-size: 12px;
                    color: var(--muted);
                    line-height: 1.4;
                }
                #controls {
                    position: absolute;
                    bottom: 16px;
                    right: 16px;
                    background: rgba(255, 255, 255, 0.88);
                    border: 1px solid var(--border);
                    border-radius: 999px;
                    box-shadow: var(--shadow);
                    backdrop-filter: blur(10px);
                    padding: 6px;
                    display: flex;
                    gap: 6px;
                    align-items: center;
                    z-index: 1000;
                }
                .control-button {
                    border: 1px solid var(--border);
                    background: linear-gradient(180deg, #ffffff 0%, #f0f4ff 100%);
                    color: var(--text);
                    border-radius: 999px;
                    padding: 8px 12px;
                    font-size: 12px;
                    font-weight: 600;
                    display: inline-flex;
                    align-items: center;
                    gap: 6px;
                    cursor: pointer;
                    transition: transform 120ms ease, box-shadow 120ms ease, border-color 120ms ease;
                }
                .control-button:hover {
                    transform: translateY(-1px);
                    box-shadow: 0 6px 16px rgba(31, 42, 61, 0.12);
                    border-color: var(--accent);
                }
                .control-button:active {
                    transform: translateY(0);
                    box-shadow: none;
                }
                .badge {
                    width: 10px;
                    height: 10px;
                    border-radius: 50%;
                    display: inline-block;
                }
                #selection-info {
                    position: absolute;
                    top: 16px;
                    left: 16px;
                    background: rgba(255, 255, 255, 0.95);
                    border: 1px solid var(--border);
                    color: var(--text);
                    padding: 12px 14px;
                    border-radius: var(--radius);
                    font-size: 12px;
                    z-index: 1000;
                    display: none;
                    box-shadow: var(--shadow);
                    max-width: 260px;
                }
                #hover-tooltip {
                    position: fixed;
                    background: rgba(17, 17, 17, 0.95);
                    color: #ffffff;
                    padding: 8px 12px;
                    border-radius: 6px;
                    font-size: 13px;
                    font-weight: 500;
                    z-index: 2000;
                    display: none;
                    pointer-events: none;
                    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.3);
                    max-width: 200px;
                    white-space: nowrap;
                }
                #hint {
                    position: absolute;
                    top: 16px;
                    left: 50%;
                    transform: translateX(-50%);
                    background: rgba(31, 42, 61, 0.85);
                    color: #ffffff;
                    padding: 8px 14px;
                    border-radius: 999px;
                    font-size: 12px;
                    z-index: 1000;
                    box-shadow: 0 10px 30px rgba(0, 0, 0, 0.25);
                    animation: fadeOut 5s ease forwards;
                }
                @keyframes fadeOut {
                    0%, 70% { opacity: 1; }
                    100% { opacity: 0; transform: translateX(-50%) translateY(-4px); }
                }
            </style>
        </head>
        <body>
            <div id="viewer"></div>

            <div id="info-panel">
                <h3><span class="dot"></span>Alprotein Structure</h3>
                <div style="font-size: 12px; margin-bottom: 6px; font-weight: 600; color: var(--text);">__STRUCT_NAME__</div>
                <div class="stat-row">
                    <span class="stat">🧬 __PIGMENT_COUNT__ pigments</span>
                    <span class="stat">🧩 __CHAIN_COUNT__ chains</span>
                </div>
                <div class="tip">Click atoms to inspect · Drag to rotate · Scroll to zoom</div>
            </div>

            <div id="controls">
                <button class="control-button" onclick="resetView()">🔄 Reset</button>
                <button class="control-button" onclick="toggleProtein()">👁️ Protein</button>
                <button class="control-button" onclick="togglePigments()">🧬 Pigments</button>
                <button class="control-button" onclick="focusOnPigments()">🎯 Focus</button>
            </div>

            <div id="selection-info">
                <div id="selection-content"></div>
            </div>

            <div id="hover-tooltip"></div>

            <div id="hint">Drag to rotate · Scroll to zoom · Shift+drag to pan</div>

            <script>
                // Initialize viewer
                let viewer = $3Dmol.createViewer('viewer', {
                    defaultcolors: $3Dmol.elementColors.rasmol,
                    backgroundColor: '#1f2a3d'
                });

                // Pigment information
                const pigmentInfo = __PIGMENT_INFO__;

                // Precompute pigment residues
                const pigmentResidues = Object.keys(pigmentInfo).map(id => parseInt(id.split('_')[2]));

                // State variables
                let showProtein = true;
                let showPigments = true;

                // Load structure
                const cifData = '__CIF_CONTENT__';
                viewer.addModel(cifData, 'cif');

                // Initial styling
                applyInitialStyling();

                let currentLabel = null;

                // Click handler - click on pigment to show label
                viewer.setClickable({}, true);
                viewer.addListener('click', function(atom) {
                    if (!atom) return;

                    // Remove previous label if it exists
                    if (currentLabel) {
                        viewer.removeLabel(currentLabel);
                        currentLabel = null;
                    }

                    // Check if clicked atom is part of a pigment
                    const clickedResi = atom.resi;
                    const isPigment = pigmentResidues.includes(clickedResi);

                    if (isPigment) {
                        // Find the pigment ID
                        const pigmentKey = Object.keys(pigmentInfo).find(key => key.includes(`_${clickedResi}_`));

                        if (pigmentKey) {
                            const pigment = pigmentInfo[pigmentKey];
                            const labelText = '🧬 ' + pigmentKey + (pigment.type ? ' (' + pigment.type + ')' : '');

                            // Add new label at the clicked atom position
                            currentLabel = viewer.addLabel(labelText, {
                                position: {x: atom.x, y: atom.y, z: atom.z},
                                backgroundColor: 'black',
                                backgroundOpacity: 0.8,
                                fontColor: 'white',
                                fontSize: 14,
                                showBackground: true,
                                alignment: 'topLeft'
                            });

                            viewer.render();
                            console.log('Label added for pigment:', pigmentKey);
                        }

                        // Also show atom info panel
                        showAtomInfo(atom);
                    } else {
                        // Still show atom info for non-pigment atoms
                        showAtomInfo(atom);
                    }
                });

                // Initial camera
                focusOnPigments();

                function applyInitialStyling() {
                    // Protein (non-pigment residues)
                    viewer.setStyle({not: {resn: ['CLA', 'CHL', 'BCL', 'CHD']}}, {
                        cartoon: {color: '#c7cedd', opacity: 0.9}
                    });

                    // Pigments with a curated palette
                    const pigmentColors = ['#2563eb', '#10b981', '#f59e0b', '#ec4899', '#8b5cf6', '#0ea5e9', '#ef4444', '#22c55e'];
                    let colorIndex = 0;

                    Object.keys(pigmentInfo).forEach(pigmentId => {
                        const color = pigmentColors[colorIndex % pigmentColors.length];
                        viewer.setStyle({resi: parseInt(pigmentId.split('_')[2])}, {
                            stick: {radius: 0.22, color: color},
                            sphere: {radius: 0.28, opacity: 0.7, color: color}
                        });
                        colorIndex++;
                    });

                    // Water and ions kept subtle
                    viewer.setStyle({resn: 'HOH'}, {sphere: {radius: 0.35, color: '#60a5fa', opacity: 0.25}});
                    viewer.setStyle({resn: ['NA', 'CL']}, {sphere: {radius: 0.7, color: '#f97316', opacity: 0.6}});

                    viewer.render();
                }

                function showAtomInfo(atom) {
                    const info = document.getElementById('selection-info');
                    const content = document.getElementById('selection-content');

                    let infoText = `<strong>Atom:</strong> ${atom.atom}<br>`;
                    infoText += `<strong>Residue:</strong> ${atom.resn} ${atom.resi}<br>`;
                    infoText += `<strong>Chain:</strong> ${atom.chain}<br>`;

                    const pigmentKey = Object.keys(pigmentInfo).find(key => key.includes(`_${atom.resi}_`));

                    if (pigmentKey) {
                        const pigment = pigmentInfo[pigmentKey];
                        infoText += `<strong>Pigment:</strong> ${pigmentKey}<br>`;
                        infoText += `<strong>Type:</strong> ${pigment.type}<br>`;
                        if (pigment.site_energy) {
                            infoText += `<strong>Site Energy:</strong> ${pigment.site_energy.toFixed(1)} cm⁻¹<br>`;
                        }
                    }

                    content.innerHTML = infoText;
                    info.style.display = 'block';

                    setTimeout(() => { info.style.display = 'none'; }, 5000);
                }

                // Control functions
                function resetView() {
                    viewer.zoomTo();
                    viewer.render();
                }

                function toggleProtein() {
                    showProtein = !showProtein;
                    if (showProtein) {
                        viewer.setStyle({not: {resn: ['CLA', 'CHL', 'BCL', 'CHD']}}, {
                            cartoon: {color: '#c7cedd', opacity: 0.9}
                        });
                    } else {
                        viewer.setStyle({not: {resn: ['CLA', 'CHL', 'BCL', 'CHD']}}, {});
                    }
                    viewer.render();
                }

                function togglePigments() {
                    showPigments = !showPigments;
                    if (showPigments) {
                        applyPigmentStyles();
                    } else {
                        Object.keys(pigmentInfo).forEach(pigmentId => {
                            viewer.setStyle({resi: parseInt(pigmentId.split('_')[2])}, {});
                        });
                    }
                    viewer.render();
                }

                function applyPigmentStyles() {
                    const pigmentColors = ['#2563eb', '#10b981', '#f59e0b', '#ec4899', '#8b5cf6', '#0ea5e9', '#ef4444', '#22c55e'];
                    let colorIndex = 0;
                    Object.keys(pigmentInfo).forEach(pigmentId => {
                        const color = pigmentColors[colorIndex % pigmentColors.length];
                        viewer.setStyle({resi: parseInt(pigmentId.split('_')[2])}, {
                            stick: {radius: 0.22, color: color},
                            sphere: {radius: 0.28, opacity: 0.7, color: color}
                        });
                        colorIndex++;
                    });
                }

                function focusOnPigments() {
                    if (pigmentResidues.length > 0) {
                        viewer.zoomTo({resi: pigmentResidues});
                    } else {
                        viewer.zoomTo();
                    }
                    viewer.render();
                }

                function hideHoverTooltip() {
                    const tooltip = document.getElementById('hover-tooltip');
                    tooltip.style.display = 'none';
                }

                // Expose viewer globally for external control
                window.alproteinViewer = viewer;
                window.alproteinPigmentInfo = pigmentInfo;
            </script>
        </body>
        </html>
        """

        html_template = html_template.replace("__WIDTH__", str(width))
        html_template = html_template.replace("__HEIGHT__", str(height))
        html_template = html_template.replace("__STRUCT_NAME__", self.structure.name)
        html_template = html_template.replace("__PIGMENT_COUNT__", str(len(self.pigment_system.pigments)))
        html_template = html_template.replace("__CHAIN_COUNT__", str(len(self.pigment_system.get_chains())))
        html_template = html_template.replace("__PIGMENT_INFO__", pigment_info)
        html_template = html_template.replace("__CIF_CONTENT__", cif_content_escaped)
        
        return html_template
    
    def generate_pdb_content(self) -> str:
        """
        Generate PDB content from the structure
        
        Returns:
            PDB format string
        """
        lines = []
        
        # Add header
        lines.append(f"HEADER    ALPROTEIN STRUCTURE                        {self.structure.name}")
        lines.append(f"TITLE     STRUCTURE FOR SITE ENERGY ANALYSIS")
        
        # Add atoms
        atom_counter = 1
        for model in self.structure.pdb:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        # Format PDB line
                        record_type = "HETATM" if residue.id[0] != " " else "ATOM"
                        
                        occupancy = getattr(atom, 'occupancy', 1.0)
                        if occupancy is None: occupancy = 1.0
                        
                        bfactor = getattr(atom, 'bfactor', 0.0)
                        if bfactor is None: bfactor = 0.0
                        
                        line = (
                            f"{record_type:<6}"
                            f"{atom_counter:>5} "
                            f"{atom.name:<4}"
                            f" "  # altloc
                            f"{residue.resname:<3} "
                            f"{chain.id}"
                            f"{residue.id[1]:>4}"
                            f" "  # icode
                            f"   "  # spaces
                            f"{atom.coord[0]:>8.3f}"
                            f"{atom.coord[1]:>8.3f}"
                            f"{atom.coord[2]:>8.3f}"
                            f"{occupancy:>6.2f}"
                            f"{bfactor:>6.2f}"
                            f"          "
                            f"{atom.element:>2}"
                        )
                        
                        lines.append(line)
                        atom_counter += 1
        
        lines.append("END")
        return "\n".join(lines)
    
    def generate_cif_content(self) -> str:
        """
        Generate mmCIF content from the structure using Bio.PDB.MMCIFIO.
        This handles large coordinates and atom counts better than PDB format.
        
        Returns:
            mmCIF format string
        """
        io_obj = MMCIFIO()
        io_obj.set_structure(self.structure.pdb)
        
        f = io.StringIO()
        io_obj.save(f)
        return f.getvalue()
    
    def get_pigment_info_js(self) -> str:
        """
        Generate JavaScript object with pigment information
        
        Returns:
            JavaScript object string
        """
        pigment_info = {}
        
        for pigment_id, pigment in self.pigment_system.pigments.items():
            pigment_info[pigment_id] = {
                'type': pigment.__class__.__name__,
                'resname': pigment.get_resname(),
                'position': pigment.get_coord().tolist(),
                'vacuum_energy': pigment.get_vacuum_energy(),
                'site_energy': None  # Will be filled when energies are calculated
            }
        
        import json
        return json.dumps(pigment_info)
    
    def generate_energy_colors(self, site_energies: Dict[str, float]) -> Dict[str, str]:
        """
        Generate color mapping based on site energies
        
        Args:
            site_energies: Dictionary mapping pigment IDs to site energies
            
        Returns:
            Dictionary mapping pigment IDs to hex color codes
        """
        if not site_energies:
            return {}
        
        # Get energy range
        energies = list(site_energies.values())
        min_energy = min(energies)
        max_energy = max(energies)
        energy_range = max_energy - min_energy
        
        # Generate colors (blue for low energy, red for high energy)
        color_map = {}
        for pigment_id, energy in site_energies.items():
            if energy_range > 0:
                # Normalize energy to 0-1 range
                normalized = (energy - min_energy) / energy_range
                
                # Convert to RGB (blue to red gradient)
                red = int(255 * normalized)
                blue = int(255 * (1 - normalized))
                green = int(128 * (1 - abs(normalized - 0.5) * 2))  # Peak at middle
                
                color_hex = f"#{red:02x}{green:02x}{blue:02x}"
            else:
                color_hex = "#888888"  # Gray if all energies are the same
            
            color_map[pigment_id] = color_hex
        
        return color_map
    
    def generate_contribution_colors(self, contributions: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
        """
        Generate color mapping based on energy contributions

        Args:
            contributions: Nested dictionary with contribution data

        Returns:
            Dictionary mapping pigment IDs to hex color codes
        """
        color_map = {}

        # Extract total shifts for coloring
        total_shifts = {}
        for pigment_id, contribs in contributions.items():
            if 'total_contribution' in contribs:
                total_shifts[pigment_id] = contribs['total_contribution']
            else:
                # Calculate total from individual contributions
                total = sum(v for k, v in contribs.items() if k != 'vacuum' and isinstance(v, (int, float)))
                total_shifts[pigment_id] = total

        return self.generate_energy_colors(total_shifts)

    def generate_domain_colors(self, domains_dict: Dict[int, List[int]]) -> Dict[str, str]:
        """
        Generate color mapping based on domain assignments

        Args:
            domains_dict: Dictionary mapping domain_id to list of pigment indices

        Returns:
            Dictionary mapping pigment IDs to hex color codes
        """
        # Define a distinct color palette for domains (up to 12 distinct colors)
        domain_palette = [
            '#2563eb',  # Blue
            '#10b981',  # Emerald
            '#f59e0b',  # Amber
            '#ec4899',  # Pink
            '#8b5cf6',  # Violet
            '#0ea5e9',  # Sky
            '#ef4444',  # Red
            '#22c55e',  # Green
            '#f97316',  # Orange
            '#06b6d4',  # Cyan
            '#a855f7',  # Purple
            '#84cc16',  # Lime
        ]

        color_map = {}

        # Get list of pigment IDs in order
        pigment_ids = list(self.pigment_system.pigments.keys())

        # Assign colors based on domain membership
        for domain_id, pigment_indices in domains_dict.items():
            # Choose color for this domain (cycle through palette if needed)
            domain_color = domain_palette[domain_id % len(domain_palette)]

            # Assign this color to all pigments in the domain
            for pig_idx in pigment_indices:
                if pig_idx < len(pigment_ids):
                    pigment_id = pigment_ids[pig_idx]
                    color_map[pigment_id] = domain_color

        return color_map
    
    def create_temp_file(self, content: str, suffix: str = ".pdb") -> str:
        """
        Create a temporary file with given content
        
        Args:
            content: File content
            suffix: File extension
            
        Returns:
            Path to temporary file
        """
        temp_fd, temp_path = tempfile.mkstemp(suffix=suffix)
        with os.fdopen(temp_fd, 'w') as temp_file:
            temp_file.write(content)
        
        self.temp_files.append(temp_path)
        return temp_path
    
    def cleanup(self):
        """Clean up temporary files"""
        for temp_file in self.temp_files:
            try:
                os.unlink(temp_file)
            except:
                pass
        self.temp_files.clear()
    
    def generate_pymol_script(self, output_path: str = None) -> str:
        """
        Generate PyMOL script for advanced visualization
        
        Args:
            output_path: Optional path to save the script
            
        Returns:
            PyMOL script string
        """
        script_lines = [
            "# Alprotein PyMOL Visualization Script",
            "# Generated automatically",
            "",
            "# Load structure",
            "load structure.pdb",
            "",
            "# Basic styling",
            "hide everything",
            "show cartoon, all",
            "color gray80, all",
            "",
            "# Style pigments",
        ]
        
        # Add pigment-specific styling
        for i, (pigment_id, pigment) in enumerate(self.pigment_system.pigments.items()):
            res_id = pigment.residue.id[1]
            chain_id = pigment.residue.get_parent().id
            
            color_name = f"pigment_color_{i}"
            script_lines.extend([
                f"set_color {color_name}, [0.{(i*37)%10}, 0.{(i*73)%10}, 0.{(i*89)%10}]",
                f"show sticks, chain {chain_id} and resi {res_id}",
                f"color {color_name}, chain {chain_id} and resi {res_id}",
            ])
        
        script_lines.extend([
            "",
            "# Final setup",
            "zoom all",
            "ray 1200, 900",
        ])
        
        script_content = "\n".join(script_lines)
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(script_content)
        
        return script_content
    
    def export_visualization_data(self, output_dir: str, site_energies: Dict[str, float] = None):
        """
        Export visualization data for external use
        
        Args:
            output_dir: Output directory
            site_energies: Optional site energy data
        """
        import os
        import json
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Export PDB file
        pdb_content = self.generate_pdb_content()
        pdb_path = os.path.join(output_dir, "structure.pdb")
        with open(pdb_path, 'w') as f:
            f.write(pdb_content)
        
        # Export pigment information
        pigment_data = {}
        for pigment_id, pigment in self.pigment_system.pigments.items():
            pigment_data[pigment_id] = {
                'type': pigment.__class__.__name__,
                'resname': pigment.get_resname(),
                'position': pigment.get_coord().tolist(),
                'vacuum_energy': pigment.get_vacuum_energy(),
                'residue_id': pigment.residue.id[1],
                'chain_id': pigment.residue.get_parent().id
            }
            
            if site_energies and pigment_id in site_energies:
                pigment_data[pigment_id]['site_energy'] = site_energies[pigment_id]
                pigment_data[pigment_id]['energy_shift'] = site_energies[pigment_id] - pigment.get_vacuum_energy()
        
        json_path = os.path.join(output_dir, "pigment_data.json")
        with open(json_path, 'w') as f:
            json.dump(pigment_data, f, indent=2)
        
        # Export PyMOL script
        pymol_script = self.generate_pymol_script()
        pymol_path = os.path.join(output_dir, "visualization.pml")
        with open(pymol_path, 'w') as f:
            f.write(pymol_script)
        
        # Export HTML viewer
        html_content = self.generate_html_view()
        html_path = os.path.join(output_dir, "viewer.html")
        with open(html_path, 'w') as f:
            f.write(html_content)
        
        print(f"Visualization data exported to: {output_dir}")
        return {
            'pdb_file': pdb_path,
            'pigment_data': json_path,
            'pymol_script': pymol_path,
            'html_viewer': html_path
        }