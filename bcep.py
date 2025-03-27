import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run Discotope3 with specified PDB file and output directory.")
    parser.add_argument("--pdb_file", type=str, required=True, help="Path to PDB file")
    parser.add_argument("--out_dir", type=str, default="test_data/test_output", help="Output directory")
    args = parser.parse_args()

    script_path = "src/discotope3_web/discotope3/main.py"
    models_dir = "src/discotope3_web/models"

    try:
        subprocess.run([
            "python3", script_path,
            "--pdb_or_zip_file", args.pdb_file,
            "--out_dir", args.out_dir,
            "--models_dir", models_dir
        ], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Discotope3: {e}")

if __name__ == "__main__":
    main()
