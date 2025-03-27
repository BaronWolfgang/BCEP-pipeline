import subprocess
import argparse
import sys

# Function to run Discotope3
def run_discotope3(pdb_file, out_dir):
    script_path = "src/discotope3_web/discotope3/main.py"
    models_dir = "src/discotope3_web/models"

    try:
        subprocess.run([
            "python3", script_path,
            "--pdb_or_zip_file", pdb_file,
            "--out_dir", out_dir,
            "--models_dir", models_dir
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Discotope3: {e}")

# Function to run Bepipred3
def run_bepipred3(fasta_file, out_dir, pred_model):
    script_path = "src/BepiPred3.0-Predictor/bepipred3_CLI.py"

    try:
        subprocess.run([
            "python3", script_path,
            "-i", fasta_file,
            "-o", out_dir,
            "-pred", pred_model
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Bepipred3: {e}")

def main():
    parser = argparse.ArgumentParser(description="Run B-cell epitope prediction tools.")
    
    # Add the tool selection argument to accept multiple tools
    parser.add_argument(
        "--tool", 
        choices=["discotope3", "bepipred3"],
        required=True,
        nargs='+',  # Allow multiple tools to be selected
        help="Select one or more tools to run (discotope3, bepipred3)"
    )

    # Common argument for output directory
    parser.add_argument(
        "--out_dir", 
        type=str, 
        default="test_data/test_output", 
        help="Output directory (default: 'test_data/test_output')"
    )
    
    # Arguments for Discotope3
    parser.add_argument(
        "--pdb", 
        type=str, 
        help="PDB file for Discotope3 (required if using discotope3)"
    )

    # Arguments for Bepipred3
    parser.add_argument(
        "--fasta", 
        type=str, 
        help="FASTA file for Bepipred3 (required if using bepipred3)"
    )

    parser.add_argument(
        "--pred", 
        type=str, 
        default="vt_pred", 
        help="Prediction model for Bepipred3 (default: 'vt_pred')"
    )

    args = parser.parse_args()

    # Loop through selected tools and run the corresponding function
    for tool in args.tool:
        if tool == "discotope3":
            if not args.pdb:
                print("Error: --pdb is required for Discotope3.")
                sys.exit(1)
            print(f"Running {tool}...")
            run_discotope3(args.pdb, args.out_dir)

        elif tool == "bepipred3":
            if not args.fasta:
                print("Error: --fasta is required for Bepipred3.")
                sys.exit(1)
            print(f"Running {tool}...")
            run_bepipred3(args.fasta, args.out_dir, args.pred)

if __name__ == "__main__":
    main()
