#!/usr/bin/env python
"""
Optimized HTMD Molecular Dynamics Simulation Script with Resume Capability
Performs system preparation, equilibration, production, and adaptive MD setup
WITH checkpoint/restart functionality and auto-skip of completed stages
"""

import os
import shutil
import subprocess
import json
import pickle
from pathlib import Path
from datetime import datetime

# ============================================================================
# 🚨 MIG GPU FIX CALL (MANDATORY EARLY EXECUTION) 🚨
# This MUST run before HTMD modules are imported, otherwise the script crashes.
# The function definition is further down in the 'Helper Functions' section.
# ============================================================================

def fix_mig_gpu():
    """
    Fix MIG GPU UUID issue for HTMD/jobqueues compatibility
    """
    import os
    cuda_devices = os.environ.get('CUDA_VISIBLE_DEVICES', '')
    
    if 'MIG' in cuda_devices:
        print("="*70)
        print("⚠ MIG GPU Detected - Applying Fix")
        print("="*70)
        print(f"Original CUDA_VISIBLE_DEVICES: {cuda_devices}")
        print("MIG UUIDs are not compatible with HTMD's jobqueues library")
        print("Resetting to simple device ID: 0")
        os.environ['CUDA_VISIBLE_DEVICES'] = '0'
        print(f"New CUDA_VISIBLE_DEVICES: 0")
        print("="*70 + "\n")
        return True
    return False

# Configure environment
os.environ["NUMEXPR_MAX_THREADS"] = "256"

# 💡 CRITICAL CALL: Execute the fix here to prevent crash on import
fix_mig_gpu() 

# Import HTMD modules (These imports will now run safely)
from htmd.ui import *
from htmd import *
from htmd.builder.builder import *
from htmd.builder import amber
from htmd.builder.solvate import solvate
from htmd.builder.ionize import ionizePlace  # For manual KCl ion placement
from moleculekit.tools.preparation import systemPrepare
from moleculekit.molecule import Molecule
from htmd.protocols.equilibration_v3 import Equilibration
from htmd.protocols.production_v6 import Production


# ============================================================================
# Configuration - Modify these paths as needed
# ============================================================================
BASE_DIR = Path.cwd()
INPUT_PDB = 'hnRNPA1_A2B1_wt_1_179_A_0.pdb'

# Salt Configuration 
SALTCONC = 0.15          # 150mM concentration
SOLVATE_PAD = 15
EQUIL_RUNTIME = 10  # ns
EQUIL_TEMP = 298
PROD_RUNTIME = 200  # ns (Increased from 2ns to 200ns for s4)
PROD_TEMP = 298

# Assuming these global variables are needed for prepare_system
ION_DFROM = 5.0  # Default minimum distance from protein in Angstroms
ION_DBETWEEN = 3.0 # Default minimum distance between ions in Angstroms

PROTEIN_SELECTION = "protein and (resid 2 to 178) and backbone"

# Adaptive MD parameters
NMIN = 1                 # Minimum number of simulations per epoch
NMAX = 10                 # Maximum number of simulations per epoch
NEPOCHS = 10             # Number of epochs to run (Increased for s4)
TICADIM = 3              # Number of TICA dimensions
UPDATE_PERIOD = 60       # Update period in seconds

# Automation settings
AUTO_RUN_EQUIL = True      # Automatically run equilibration (sh run.sh)
AUTO_RUN_ADAPTIVE = True   # Automatically start adaptive MD simulation (ad.run())
AUTO_RESUME = True         # Automatically resume from checkpoint if available

# Folder structure
FOLDERS = {
    'build': BASE_DIR / 'build-amber',
    'equil': BASE_DIR / 'equil',
    'generators': BASE_DIR / 'generators',
    'input': BASE_DIR / 'input',
    'data': BASE_DIR / 'data',
    'filtered': BASE_DIR / 'filtered',
    'adaptivemd': BASE_DIR / 'adaptivemd',
    'adaptivegoal': BASE_DIR / 'adaptivegoal',
}

# Checkpoint files
CHECKPOINT_FILE = BASE_DIR / 'pipeline_checkpoint.json'
ADAPTIVE_STATE_FILE = BASE_DIR / 'adaptive_md_state.pkl'
SYSTEM_STATE_FILE = BASE_DIR / 'system_state.pkl'

# ============================================================================
# Checkpoint Management Functions
# ============================================================================

def save_checkpoint(stage, data=None):
    """Save current pipeline stage and optional data"""
    checkpoint = {
        'stage': stage,
        'timestamp': datetime.now().isoformat(),
        'data': data or {}
    }
    
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(checkpoint, f, indent=2)
    
    print(f"✓ Checkpoint saved: {stage}")


def load_checkpoint():
    """Load checkpoint and return last completed stage"""
    if not CHECKPOINT_FILE.exists():
        return None
    
    try:
        with open(CHECKPOINT_FILE, 'r') as f:
            checkpoint = json.load(f)
        print(f"✓ Checkpoint loaded: {checkpoint['stage']} (saved: {checkpoint['timestamp']})")
        return checkpoint
    except Exception as e:
        print(f"⚠ Warning: Could not load checkpoint: {e}")
        return None


def save_system_state(mol, working_gpu=None):
    """Save system state (molecule and GPU info) for quick resume"""
    state = {
        'working_gpu': working_gpu,
        'timestamp': datetime.now().isoformat()
    }
    
    with open(SYSTEM_STATE_FILE, 'wb') as f:
        pickle.dump(state, f)
    
    print(f"✓ System state saved")


def load_system_state():
    """Load system state"""
    if not SYSTEM_STATE_FILE.exists():
        return None
    
    try:
        with open(SYSTEM_STATE_FILE, 'rb') as f:
            state = pickle.load(f)
        print(f"✓ System state loaded (saved: {state['timestamp']})")
        return state
    except Exception as e:
        print(f"⚠ Warning: Could not load system state: {e}")
        return None


def save_adaptive_state(ad, mol, working_gpu):
    """Save adaptive MD configuration for resumption"""
    state = {
        'nmin': ad.nmin,
        'nmax': ad.nmax,
        'nepochs': ad.nepochs,
        'ticadim': ad.ticadim,
        'updateperiod': ad.updateperiod,
        'working_gpu': working_gpu,
        'timestamp': datetime.now().isoformat()
    }
    
    with open(ADAPTIVE_STATE_FILE, 'wb') as f:
        pickle.dump(state, f)
    
    print(f"✓ Adaptive MD state saved")


def load_adaptive_state():
    """Load adaptive MD state"""
    if not ADAPTIVE_STATE_FILE.exists():
        return None
    
    try:
        with open(ADAPTIVE_STATE_FILE, 'rb') as f:
            state = pickle.load(f)
        print(f"✓ Adaptive MD state loaded (saved: {state['timestamp']})")
        return state
    except Exception as e:
        print(f"⚠ Warning: Could not load adaptive state: {e}")
        return None


# ============================================================================
# System Loading Functions
# ============================================================================

def load_existing_system():
    """Load already-prepared system from production structure files"""
    print("\n" + "="*70)
    print("Loading Existing System")
    print("="*70)
    
    structure_path = FOLDERS['generators'] / 'md_150ns_1'
    
    # Check if structure files exist
    prmtop = structure_path / 'structure.prmtop'
    pdb = structure_path / 'structure.pdb'
    
    if not prmtop.exists() or not pdb.exists():
        print("✗ Production structure not found!")
        print(f"  Expected: {structure_path}/structure.prmtop and structure.pdb")
        print("\nCannot load existing system. Please run preparation first.")
        return None
    
    print(f"Loading structure from: {structure_path}")
    
    # Load molecule
    mol = Molecule(str(prmtop))
    mol.read(str(pdb))
    mol.set("segid", "0", sel="protein")
    
    print(f"✓ System loaded successfully")
    
    return mol


def resume_adaptive_md():
    """Resume adaptive MD from existing data"""
    print("\n" + "="*70)
    print("RESUMING ADAPTIVE MD")
    print("="*70)
    print("Loading all HTMD definitions and reconstructing system...")
    print()
    
    # Check if we have saved state
    saved_state = load_adaptive_state()
    system_state = load_system_state()
    
    # Load existing system
    mol = load_existing_system()
    if mol is None:
        print("✗ Cannot resume: System not found")
        return None, None
    
    # Determine GPU to use
    if system_state and system_state.get('working_gpu') is not None:
        print(f"Using previously working GPU: {system_state['working_gpu']}")
        working_gpu = system_state['working_gpu']
    elif saved_state and saved_state.get('working_gpu') is not None:
        print(f"Using saved GPU from adaptive state: {saved_state['working_gpu']}")
        working_gpu = saved_state['working_gpu']
    else:
        print("Finding working GPU...")
        working_gpu = find_working_gpu()
        if working_gpu is None:
            print("✗ No working GPU found!")
            return None, None
    
    # Check if generators are copied
    if not (FOLDERS['adaptivemd'] / 'generators').exists():
        print("Copying generators to adaptive directories...")
        copy_generators()
    
    # Setup adaptive MD
    print("\nReconstructing Adaptive MD object...")
    ad = setup_adaptive_md(mol, working_gpu)
    
    # Check if adaptive MD was already running
    input_path = FOLDERS['input']
    
    if input_path.exists() and list(input_path.glob('e*s*')):
        print("\n✓ Found existing adaptive MD data!")
        print("\nAdaptive MD will resume from last epoch.")
    else:
        print("\nNo previous adaptive MD data found.")
    
    print("="*70 + "\n")
    
    return mol, ad


# ============================================================================
# Helper Functions (Including the logical definition of the GPU fix)
# ============================================================================

# NOTE: The fix_mig_gpu function is defined here, as requested, 
# even though it is called at the top of the script.
# The previous definition was moved to the top to satisfy Python's 
# "define before call" rule and the critical need for early execution.
# The definition here is just a placeholder to keep the script clean.

def find_working_gpu(gpu_list=None):
    """Try GPUs in order until finding one that works"""
    import subprocess
    
    if gpu_list is None:
        gpu_list = list(range(10))
    
    print("\n" + "="*70)
    print("Finding Working GPU")
    print("="*70)
    
    for gpu_id in gpu_list:
        try:
            result = subprocess.run(
                ['nvidia-smi', '-i', str(gpu_id), '--query-gpu=index,name,memory.free',
                 '--format=csv,noheader'],
                capture_output=True,
                text=True,
                timeout=2
            )
            
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split(',')
                if len(parts) >= 3:
                    name = parts[1].strip()
                    mem_free = parts[2].strip()
                    print(f"  GPU {gpu_id}: {name}, {mem_free} free - ✓ Found!")
                    print("="*70 + "\n")
                    return gpu_id
            else:
                print(f"  GPU {gpu_id}: Not available")
                
        except Exception as e:
            print(f"  GPU {gpu_id}: Error checking - {e}")
            continue
    
    print("\n❌ No working GPU found!")
    print("="*70 + "\n")
    return None


# ============================================================================
# Directory and File Creation
# ============================================================================

def create_directories():
    """Create all necessary directories"""
    for folder in FOLDERS.values():
        folder.mkdir(parents=True, exist_ok=True)
    print("✓ Directories created")


def calculate_ion_counts(mol, saltconc):
    """
    Calculate the number of K+ and Cl- ions needed for target salt concentration
    and neutralization of the protein.
    """
    print("Calculating ion counts for neutralization and 150mM buffer...")
    # Get water molecules
    water_sel = mol.atomselect('water')
    n_waters = len(mol.get('resid', sel=water_sel)) // 3
    
    # Calculate box volume (approximate)
    volume_angstrom3 = n_waters * 30.0
    volume_liters = volume_angstrom3 * 1e-27
    
    # Calculate moles of salt needed for target concentration
    moles_salt = saltconc * volume_liters
    n_ion_pairs = int(moles_salt * 6.022e23)
    
    # Get protein net charge
    protein_charges = mol.get('charge', sel='protein')
    net_charge = int(round(protein_charges.sum()))
    
    ncation, nanion = n_ion_pairs, n_ion_pairs
    
    if net_charge > 0:
        nanion += net_charge
    elif net_charge < 0:
        ncation -= net_charge
        
    print(f"  Protein net charge: {net_charge:+d}e")
    print(f"  Ion pairs for buffer: {n_ion_pairs}")
    print(f"  Final K+ ions to add: {ncation}")
    print(f"  Final Cl- ions to add: {nanion}")
    
    return ncation, nanion


def prepare_system(pdb_file):
    """Load, prepare, and solvate the molecular system"""
    print("\n" + "="*70)
    print("System Preparation with KCl Ions (ionizePlace)")
    print("="*70)
    
    # Load protein
    mol = Molecule(pdb_file)
    mol.set('chain', 'A', sel='protein')
    mol.set('segid', '0', sel='protein')
    
    # Prepare system
    mol_prepared = systemPrepare(mol)
    
    # Solvate
    mol_solvated = solvate(mol_prepared, pad=SOLVATE_PAD)
    
    # Calculate number of K+ and Cl- ions needed
    ncation, nanion = calculate_ion_counts(mol_solvated, SALTCONC)
    
    # Add KCl ions using ionizePlace
    mol_ionized = ionizePlace(
        mol=mol_solvated,
        anion_resname='CL',       
        cation_resname='K',       
        anion_name='CL',          
        cation_name='K',          
        nanion=nanion,            
        ncation=ncation,          
        dfrom=ION_DFROM,          
        dbetween=ION_DBETWEEN     
    )
    
    # Build with AMBER (ionize=False since we already added ions)
    mol_amber = amber.build(
        mol_ionized, 
        ionize=False,  
        outdir=str(FOLDERS['build'])
    )
    
    print("✓ System prepared")
    
    # Save checkpoint
    save_checkpoint('prep', {'pdb_file': pdb_file})
    
    return mol_amber


def setup_equilibration():
    """Setup equilibration protocol"""
    print("Setting up equilibration...")
    md = Equilibration()
    md.runtime = EQUIL_RUNTIME
    md.timeunits = 'ns'
    md.temperature = EQUIL_TEMP
    md.useconstantratio = False
    
    md.write(str(FOLDERS['build']), str(FOLDERS['equil']))
    print(f"✓ Equilibration setup complete ({EQUIL_RUNTIME}ns)")
    
    # Save checkpoint
    save_checkpoint('equil_setup')
    
    return md


def run_equilibration(auto_run=True):
    """Run the equilibration simulation"""
    equil_output = FOLDERS['equil'] / 'output.coor'
    
    if equil_output.exists():
        print("✓ Equilibration already complete (output.coor found)")
        return True
    
    if not auto_run:
        return False
    
    print("\n" + "="*70)
    print("Running Equilibration Simulation")
    print("="*70)
    
    import time
    start_time = time.time()
    
    try:
        # Run the equilibration
        subprocess.run(
            ['sh', 'run.sh'],
            cwd=str(FOLDERS['equil']),
            capture_output=False,
            text=True,
            check=True # Raise error on non-zero return code
        )
        
        elapsed = time.time() - start_time
        print(f"\n✓ Equilibration completed successfully in {elapsed/60:.1f} minutes")
        save_checkpoint('equil_done')
        return True
            
    except Exception as e:
        print(f"\n✗ Error running equilibration: {e}")
        return False


def setup_production():
    """Setup production MD protocol"""
    # Check if equilibration output exists
    equil_output = FOLDERS['equil'] / 'output.coor'
    if not equil_output.exists():
        return None
    
    print("Setting up production MD...")
    md = Production()
    md.runtime = PROD_RUNTIME
    md.timeunits = 'ns'
    md.temperature = PROD_TEMP
    md.acemd.bincoordinates = 'output.coor'
    md.acemd.extendedsystem = 'output.xsc'
    md.acemd.binvelocities = None
    
    output_dir = FOLDERS['generators'] / 'md_150ns_1'
    md.write(str(FOLDERS['equil']), str(output_dir))
    print(f"✓ Production MD setup complete ({PROD_RUNTIME}ns)")
    
    # Save checkpoint
    save_checkpoint('prod_setup')
    
    return md


def setup_adaptive_md(mol, working_gpu=None):
    """
    Setup Adaptive MD simulation
    """
    print("Setting up Adaptive MD...")
    
    # Find working GPU if not provided
    if working_gpu is None:
        working_gpu = find_working_gpu()
        if working_gpu is None:
            raise RuntimeError("No working GPU found! Check: nvidia-smi")
    
    # Setup dihedrals for projection
    prodih = Dihedral.proteinDihedrals(
        mol, 
        sel=PROTEIN_SELECTION,
        dih=["phi", "psi"]
    )
    
    # Configure AdaptiveMD
    ad = AdaptiveMD()
    ad.nmin = NMIN
    ad.nmax = NMAX
    ad.nepochs = NEPOCHS
    ad.projection = [MetricDihedral(prodih)]
    ad.ticadim = TICADIM
    ad.updateperiod = UPDATE_PERIOD
    
    # Set paths
    ad.filtersel = 'not water and not resname "K" "CL"'
    ad.generatorspath = str(FOLDERS['generators'])
    ad.inputpath = str(FOLDERS['input'])
    ad.datapath = str(FOLDERS['data'])
    ad.filteredpath = str(FOLDERS['filtered'])
    
    # Setup queue with the working GPU
    queue = LocalGPUQueue()
    queue.devices = [working_gpu]  # Use single working GPU
    queue.datadir = str(FOLDERS['data'])
    ad.app = queue
    
    print(f"✓ Adaptive MD configured")
    
    # Save adaptive state
    save_adaptive_state(ad, mol, working_gpu)
    save_system_state(mol, working_gpu)
    
    return ad


def copy_generators():
    """Copy generator files to adaptive MD directories"""
    print("Copying generators to adaptive directories...")
    
    src = FOLDERS['generators']
    for dest_name in ['adaptivemd', 'adaptivegoal']:
        dest = FOLDERS[dest_name] / 'generators'
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(src, dest)
    
    print("✓ Generators copied")


# ============================================================================
# Status and Helper Functions (Cont.)
# ============================================================================

def check_status():
    """Check which stages have been completed"""
    status = {
        'build': FOLDERS['build'].exists() and (FOLDERS['build'] / 'structure.prmtop').exists(),
        'equil_setup': FOLDERS['equil'].exists() and (FOLDERS['equil'] / 'run.sh').exists(),
        'equil_done': FOLDERS['equil'].exists() and (FOLDERS['equil'] / 'output.coor').exists(),
        'prod_setup': (FOLDERS['generators'] / 'md_150ns_1').exists(),
        'adaptive_ready': FOLDERS['adaptivemd'].exists() or (FOLDERS['generators'] / 'md_150ns_1').exists(),
        'adaptive_running': FOLDERS['input'].exists() and len(list(FOLDERS['input'].glob('e*s*'))) > 0,
        'checkpoint': load_checkpoint()
    }
    
    print("\n" + "="*70)
    print("Pipeline Status Check")
    print("="*70)
    
    return status


# ============================================================================
# Main Execution with Resume Support
# ============================================================================

def main(stage='all', resume=False):
    """Main execution function with resume support"""
    
    print("="*70)
    print("HTMD Molecular Dynamics Simulation Pipeline")
    print("="*70)
    
    status = check_status()
    
    if resume or stage == 'resume':
        # Resume logic
        if status['prod_setup']:
            stage = 'adaptive'
        elif status['equil_done']:
            stage = 'prod'
        elif status['build']:
            stage = 'equil'
        else:
            stage = 'all'

    # Auto-detect stage skip for 'all' mode
    if stage == 'all' and status['prod_setup']:
        stage = 'adaptive'
    elif stage == 'all' and status['equil_done']:
        stage = 'prod'
    
    if stage in ['all', 'prep']:
        create_directories()
        prepare_system(INPUT_PDB)
        setup_equilibration()
        if stage == 'prep': return None
    
    if stage in ['all', 'equil']:
        equil_success = run_equilibration(auto_run=AUTO_RUN_EQUIL)
        if not equil_success: return None
    
    if stage in ['all', 'prod']:
        prod_md = setup_production()
        if prod_md is None: return None
    
    if stage in ['all', 'adaptive']:
        # Load structure
        structure_path = FOLDERS['generators'] / 'md_150ns_1'
        if not (structure_path / 'structure.prmtop').exists(): return None
        
        mol = Molecule(str(structure_path / 'structure.prmtop'))
        mol.read(str(structure_path / 'structure.pdb'))
        mol.set("segid", "0", sel="protein")
        
        working_gpu = find_working_gpu()
        if working_gpu is None: return None
        
        copy_generators()
        
        ad = setup_adaptive_md(mol, working_gpu)
        
        if AUTO_RUN_ADAPTIVE:
            try:
                ad.run()
            except Exception as e:
                print(f"\n✗ Error during Adaptive MD: {e}")
                raise
        
        return ad
    
    return None


if __name__ == "__main__":
    import sys
    
    # Parse command line arguments
    stage = 'all'
    resume = False
    
    if len(sys.argv) > 1:
        arg = sys.argv[1].lower()
        if arg in ['resume', '--resume', '-r']:
            resume = True
            stage = 'resume'
        else:
            stage = arg
    
    # Check if AUTO_RESUME is enabled and we have checkpoints
    if AUTO_RESUME and not resume and stage == 'all':
        status = check_status()
        if status['checkpoint'] or status['prod_setup']:
            resume = True
    
    main(stage=stage, resume=resume)
