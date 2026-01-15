
import os
import numpy as np
import pandas as pd
from pathlib import Path
import json

# Try to import RDKit for ADMET calculations
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

class DFTDataCollector:
    """
    Collects and aggregates DFT calculation results from wavefunction files (.npz)
    into a comprehensive dataset for analysis, including ADMET properties.
    """
    
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        
    def collect_data(self):
        """
        Scans the output directory for optimized wavefunction files and compiles stats.
        Returns a Pandas DataFrame with DFT and ADMET properties.
        """
        data_records = []
        
        # Find all wavefunction files recursively
        # Pattern: */wavefunctions/*_optimized_wavefunction.npz
        wfn_files = list(self.output_dir.rglob("*_optimized_wavefunction.npz"))
        
        print(f"Found {len(wfn_files)} wavefunction files in {self.output_dir}")
        
        for wfn_file in wfn_files:
            try:
                record = self._process_single_file(wfn_file)
                if record:
                    data_records.append(record)
            except Exception as e:
                print(f"Error processing {wfn_file.name}: {e}")
                
        # Create DataFrame
        df = pd.DataFrame(data_records)
        
        # Reorder columns for logical flow if data exists
        if not df.empty:
            preferred_order = [
                # Identification
                'Molecule_Name', 'SMILES', 'Formula', 
                # DFT Electronic
                'Total_Energy_Ha', 'HOMO_eV', 'LUMO_eV', 'Gap_eV', 
                'Dipole_Moment_Debye', 'Polarizability_au3', 'Volume_A3', 'PSA_A2',
                'Hardness_eV', 'Electrophilicity_Index', 
                'Max_Fukui_f_plus', 'Max_Fukui_f_minus', 'Max_Spin_Density',
                'ESP_Min_au', 'ESP_Max_au', 'ESP_Variance',
                # ADMET Properties
                'MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotatableBonds', 
                'AromaticRings', 'FractionCsp3', 'QED',
                # Drug-likeness Filters
                'Lipinski_Pass', 'Lipinski_Violations', 
                'Veber_Pass', 'Ghose_Pass', 'Egan_Pass', 'Muegge_Pass', 'LeadLike_Pass',
                # ADMET Predictions
                'BBB_Penetration', 'hERG_Risk', 'Solubility', 'CYP_Inhibition',
                'Aggregator_Risk', 'PAINS_Alerts', 'Oral_Bioavailability'
            ]
            
            # Get columns that actually exist in the dataframe
            existing_cols = [c for c in preferred_order if c in df.columns]
            # Add remaining columns
            remaining_cols = [c for c in df.columns if c not in existing_cols]
            
            df = df[existing_cols + remaining_cols]
            
        return df
    
    def _calculate_admet_properties(self, smiles):
        """Calculate ADMET properties from SMILES string."""
        if not RDKIT_AVAILABLE or not smiles:
            return {}
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {}
            
            mol_no_h = Chem.RemoveHs(mol)
            mol_with_h = Chem.AddHs(mol)
            
            # Basic descriptors
            props = {
                'MW': Descriptors.MolWt(mol_no_h),
                'LogP': Descriptors.MolLogP(mol_no_h),
                'TPSA': Descriptors.TPSA(mol_no_h),
                'HBD': Descriptors.NumHDonors(mol_no_h),
                'HBA': Descriptors.NumHAcceptors(mol_no_h),
                'RotatableBonds': Descriptors.NumRotatableBonds(mol_no_h),
                'AromaticRings': Descriptors.NumAromaticRings(mol_no_h),
                'FractionCsp3': Descriptors.FractionCSP3(mol_no_h),
                'HeavyAtoms': mol_no_h.GetNumHeavyAtoms(),
                'RingCount': Descriptors.RingCount(mol_no_h),
                'QED': QED.qed(mol_no_h)
            }
            
            # Lipinski Rule of 5
            lipinski_violations = 0
            if props['MW'] > 500: lipinski_violations += 1
            if props['LogP'] > 5: lipinski_violations += 1
            if props['HBD'] > 5: lipinski_violations += 1
            if props['HBA'] > 10: lipinski_violations += 1
            props['Lipinski_Violations'] = lipinski_violations
            props['Lipinski_Pass'] = lipinski_violations <= 1
            
            # Veber rules
            veber_pass = props['TPSA'] <= 140 and props['RotatableBonds'] <= 10
            props['Veber_Pass'] = veber_pass
            
            # Ghose filter
            ghose_pass = (160 <= props['MW'] <= 480 and 
                         -0.4 <= props['LogP'] <= 5.6 and 
                         20 <= props['HeavyAtoms'] <= 70)
            props['Ghose_Pass'] = ghose_pass
            
            # Egan filter
            egan_pass = props['TPSA'] <= 131.6 and -1 <= props['LogP'] <= 6
            props['Egan_Pass'] = egan_pass
            
            # Muegge filter
            muegge_pass = (200 <= props['MW'] <= 600 and 
                          -2 <= props['LogP'] <= 5 and 
                          props['TPSA'] <= 150 and 
                          props['RingCount'] <= 7)
            props['Muegge_Pass'] = muegge_pass
            
            # Lead-like
            leadlike_pass = (250 <= props['MW'] <= 350 and 
                            props['LogP'] <= 4 and 
                            props['RotatableBonds'] <= 7)
            props['LeadLike_Pass'] = leadlike_pass
            
            # BBB Penetration (Egan/Clark rules)
            bbb_score = 0
            if props['TPSA'] < 90: bbb_score += 1
            if 1 <= props['LogP'] <= 3: bbb_score += 1
            if props['MW'] < 450: bbb_score += 1
            if props['HBD'] <= 2: bbb_score += 1
            props['BBB_Penetration'] = 'High' if bbb_score >= 4 else ('Moderate' if bbb_score >= 2 else 'Low')
            
            # hERG Risk (Aronov model)
            n_basic = len(mol_with_h.GetSubstructMatches(Chem.MolFromSmarts('[#7;+,H1,H2]')))
            herg_risk = 0
            if props['LogP'] > 3: herg_risk += 1
            if props['LogP'] > 4: herg_risk += 1
            if n_basic > 0: herg_risk += 1
            if n_basic > 1 and props['LogP'] > 3: herg_risk += 1
            if props['MW'] > 400 and props['TPSA'] < 75: herg_risk += 1
            props['hERG_Risk'] = 'High' if herg_risk >= 3 else ('Moderate' if herg_risk >= 1 else 'Low')
            
            # Solubility (GSK 4/400 rule)
            sol_score = 0
            if props['LogP'] < 3: sol_score += 2
            elif props['LogP'] < 4: sol_score += 1
            if props['MW'] < 400: sol_score += 1
            if props['AromaticRings'] <= 3: sol_score += 1
            if props['TPSA'] > 50: sol_score += 1
            props['Solubility'] = 'Good' if sol_score >= 4 else ('Moderate' if sol_score >= 2 else 'Poor')
            
            # CYP Inhibition risk
            cyp_risk = 0
            if props['LogP'] > 3: cyp_risk += 1
            if props['AromaticRings'] > 2: cyp_risk += 1
            if props['MW'] > 350: cyp_risk += 1
            props['CYP_Inhibition'] = 'High' if cyp_risk >= 3 else ('Moderate' if cyp_risk >= 1 else 'Low')
            
            # Aggregator Risk (Shoichet model)
            agg_risk = 0
            if props['LogP'] > 4: agg_risk += 2
            elif props['LogP'] > 3: agg_risk += 1
            if props['AromaticRings'] > 3: agg_risk += 2
            elif props['AromaticRings'] > 2: agg_risk += 1
            if props['TPSA'] < 50: agg_risk += 1
            if props['HeavyAtoms'] > 35: agg_risk += 1
            if props['FractionCsp3'] < 0.2: agg_risk += 1
            props['Aggregator_Risk'] = 'High' if agg_risk >= 4 else ('Medium' if agg_risk >= 2 else 'Low')
            
            # PAINS check (simplified)
            pains_patterns = {
                'quinone': '[#6]1([#8])=[#6][#6]([#8])=[#6][#6]=[#6]1',
                'catechol': 'c1cc(O)c(O)cc1',
                'michael_acceptor': '[CH2]=[CH][C,c](=O)',
                'nitro_aromatic': 'c[N+](=O)[O-]',
                'aldehyde': '[CH]=O',
            }
            pains_alerts = []
            for name, smarts in pains_patterns.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol_with_h.HasSubstructMatch(pattern):
                    pains_alerts.append(name)
            props['PAINS_Alerts'] = ','.join(pains_alerts) if pains_alerts else 'None'
            
            # Oral Bioavailability (combined Lipinski + Veber)
            props['Oral_Bioavailability'] = 'Yes' if (props['Lipinski_Pass'] and props['Veber_Pass']) else 'No'
            
            # Remove intermediate keys not needed in CSV
            props.pop('HeavyAtoms', None)
            props.pop('RingCount', None)
            
            return props
            
        except Exception as e:
            print(f"ADMET calculation error: {e}")
            return {}
    
    def _process_single_file(self, file_path):
        """Extracts properties from a single .npz file."""
        try:
            data = np.load(file_path, allow_pickle=True)
            
            # infer molecule name from directory or filename
            # Structure usually: .../MoleculeName/wavefunctions/MoleculeName_optimized_wfn.npz
            mol_name = file_path.parent.parent.name
            
            record = {
                'Molecule_Name': mol_name,
                'File_Path': str(file_path)
            }
            
            # --- 1. Basic Electronic Structure ---
            if 'energy' in data:
                record['Total_Energy_Ha'] = float(data['energy'])
            if 'dipole_magnitude' in data:
                record['Dipole_Moment_Debye'] = float(data['dipole_magnitude'])
            
            # HOMO/LUMO/Gap
            if 'mo_energy' in data and 'mo_occ' in data:
                mo_energy = data['mo_energy']
                mo_occ = data['mo_occ']
                
                # Identify HOMO/LUMO
                occ_idx = np.where(mo_occ > 0)[0]
                if len(occ_idx) > 0:
                    homo_idx = occ_idx[-1]
                    lumo_idx = homo_idx + 1 if homo_idx + 1 < len(mo_energy) else homo_idx
                    
                    homo_ev = mo_energy[homo_idx] * 27.2114
                    lumo_ev = mo_energy[lumo_idx] * 27.2114
                    gap_ev = lumo_ev - homo_ev
                    
                    record['HOMO_eV'] = homo_ev
                    record['LUMO_eV'] = lumo_ev
                    record['Gap_eV'] = gap_ev
            
            # --- 2. Conceptual DFT Descriptors ---
            # (Calculated fresh to ensure consistency, or read if available)
            # Keys might vary, so let's recalculate from HOMO/LUMO if missing
            if 'chemical_hardness' in data:
                record['Hardness_eV'] = float(data['chemical_hardness'])
            elif 'HOMO_eV' in record:
                ip = -record['HOMO_eV']
                ea = -record['LUMO_eV']
                record['Hardness_eV'] = (ip - ea) / 2
                
            if 'electrophilicity_index' in data:
                record['Electrophilicity_Index'] = float(data['electrophilicity_index'])
            elif 'HOMO_eV' in record:
                ip = -record['HOMO_eV']
                ea = -record['LUMO_eV']
                neg = (ip + ea) / 2
                hardness = (ip - ea) / 2
                record['Electrophilicity_Index'] = (neg**2)/(2*hardness) if hardness > 0 else 0
                
            if 'nucleophilicity_index' in data:
                record['Nucleophilicity_Index'] = float(data['nucleophilicity_index'])
                
            # --- 3. Molecular Properties ---
            if 'polarizability' in data:
                record['Polarizability_au3'] = float(data['polarizability'])
            if 'molecular_volume' in data:
                record['Volume_A3'] = float(data['molecular_volume'])
            if 'polar_surface_area' in data:
                record['PSA_A2'] = float(data['polar_surface_area'])
                
            # --- 4. Reactivity & Toxicity Indicators (Max Values) ---
            # For arrays like Fukui, we usually want the Maximum value to spot "hotspots"
            if 'fukui_plus' in data:
                fp = data['fukui_plus']
                record['Max_Fukui_f_plus'] = float(np.max(fp)) if fp.size > 0 else 0
                
            if 'fukui_minus' in data:
                fm = data['fukui_minus']
                record['Max_Fukui_f_minus'] = float(np.max(np.abs(fm))) if fm.size > 0 else 0 # abs for f- just in case
                
            if 'fukui_radical' in data:
                fr = data['fukui_radical']
                record['Max_Fukui_f_radical'] = float(np.max(fr)) if fr.size > 0 else 0
                
            if 'spin_densities' in data:
                sd = data['spin_densities']
                record['Max_Spin_Density'] = float(np.max(np.abs(sd))) if sd.size > 0 else 0
                
            # --- 5. Electrostatic Potential Stats ---
            if 'esp_min' in data:
                record['ESP_Min_au'] = float(data['esp_min'])
            if 'esp_max' in data:
                record['ESP_Max_au'] = float(data['esp_max'])
            if 'esp_variance' in data:
                record['ESP_Variance'] = float(data['esp_variance'])
            if 'esp_positive_avg' in data:
                record['ESP_Avg_Pos_au'] = float(data['esp_positive_avg'])
            if 'esp_negative_avg' in data:
                record['ESP_Avg_Neg_au'] = float(data['esp_negative_avg'])
                
            # --- 6. Get SMILES and calculate ADMET ---
            smiles = None
            log_file = file_path.parent.parent / "optimization_output.log"
            if log_file.exists():
                try:
                    with open(log_file, 'r') as f:
                        for line in f:
                            if "SMILES:" in line:
                                smiles = line.split("SMILES:")[1].strip()
                                record['SMILES'] = smiles
                            if "Formula:" in line:
                                record['Formula'] = line.split("Formula:")[1].strip()
                            if smiles:
                                break
                except:
                    pass
            
            # --- 7. Calculate ADMET Properties ---
            if smiles:
                admet_props = self._calculate_admet_properties(smiles)
                record.update(admet_props)
            
            return record
            
        except Exception as e:
            print(f"Failed to read {file_path}: {e}")
            return None

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Collect DFT results into CSV")
    parser.add_argument("--output_dir", default="dft_data", help="Root directory of DFT outputs")
    parser.add_argument("--csv_name", default="dft_analysis_results.csv", help="Output CSV filename")
    
    args = parser.parse_args()
    
    collector = DFTDataCollector(args.output_dir)
    df = collector.collect_data()
    
    if not df.empty:
        out_path = Path(args.csv_name).resolve()
        df.to_csv(out_path, index=False)
        print(f"✅ Successfully saved consolidated data to: {out_path}")
        print(f"   Total Layout: {df.shape[0]} molecules x {df.shape[1]} properties")
    else:
        print("⚠️ No data collected. Check your output directory path.")

if __name__ == "__main__":
    main()

