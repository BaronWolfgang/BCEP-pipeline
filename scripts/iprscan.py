#!/usr/bin/env python3
import requests
import time
import os
import argparse

# Allowed InterProScan applications
ALLOWED_APPL = [
    "NCBIfam", "SFLD", "Phobius", "SignalP", "SignalP_EUK",
    "SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE", "SuperFamily",
    "Panther", "Gene3d", "HAMAP", "PrositeProfiles", "PrositePatterns",
    "Coils", "SMART", "CDD", "PRINTS", "PfamA", "MobiDBLite",
    "PIRSF", "TMHMM", "AntiFam", "FunFam", "PIRSR"
]

def submit_ipr(fasta_path, email, appl=ALLOWED_APPL, stype='p', outfolder='ipr_out', poll_interval=5):
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    # Validate appl options
    for a in appl:
        if a not in ALLOWED_APPL:
            raise ValueError(f"Invalid appl: {a}. Must be one of {ALLOWED_APPL}")

    # Read FASTA sequence
    with open(fasta_path) as f:
        sequence = ''.join(line.strip() for line in f if not line.startswith('>'))

    # Submit job
    url_submit = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
    data = {
        'email': email,
        'sequence': sequence,
        'stype': stype,
        'appl': ','.join(appl),
        'goterms': 'false',
        'pathways': 'false',
        'title': 'iprscan_python_minimal'
    }

    r = requests.post(url_submit, data=data)
    if r.status_code != 200:
        raise RuntimeError(f"Submission failed: {r.text}")

    job_id = r.text.strip()
    print(f"Job submitted. Job ID: {job_id}")

    # Poll job until completion
    url_result = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/tsv"
    while True:
        r = requests.get(url_result, headers={'Accept': 'text/tab-separated-values'})
        if r.status_code == 200:
            break
        print("Job still running, waiting...")
        time.sleep(poll_interval)

    # Save result
    basename = os.path.splitext(os.path.basename(fasta_path))[0]  # remove extension
    outpath = os.path.join(outfolder, f"{basename}.tsv")
    with open(outpath, 'wb') as f:
        f.write(r.content)

    print(f"Job completed. Result saved to: {outpath}")
    return outpath

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit a FASTA file to InterProScan")
    parser.add_argument("--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("--email", required=True, help="Your email for job submission")
    parser.add_argument("--outfolder", default="temp/ipr_output", help="Folder to save results")
    parser.add_argument("--stype", default="p", choices=["p", "n"], help="Sequence type: 'p' for protein, 'n' for nucleotide")
    parser.add_argument("--appl", nargs="+", default=ALLOWED_APPL, help=f"Analyses to run. Allowed: {ALLOWED_APPL}")
    args = parser.parse_args()

    submit_ipr(
        fasta_path=args.fasta,
        email=args.email,
        outfolder=args.outfolder,
        stype=args.stype,
        appl=args.appl
    )
