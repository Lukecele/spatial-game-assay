#!/usr/bin/env python3
"""
download_data.py
================
Download and unpack authentic 10x Genomics Visium CytAssist spatial transcriptomics data
for Human Lung Squamous Cell Carcinoma.
Dataset URL: https://www.10xgenomics.com/resources/datasets/human-lung-cancer-ffpe-2-standard
"""

import os
import sys
import tarfile
import urllib.request
import time

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

URLS = {
    "spatial": (
        "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/"
        "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/"
        "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_spatial.tar.gz"
    ),
    "matrix": (
        "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/"
        "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma/"
        "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.tar.gz"
    )
}

HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
}


def download_file(url, dest_path):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req) as resp, open(dest_path, "wb") as out_file:
        total_size = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 1024  # 1 MB chunks
        start_time = time.time()
        
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            out_file.write(chunk)
            downloaded += len(chunk)
            if total_size > 0:
                percent = min(100.0, downloaded * 100.0 / total_size)
                mb_down = downloaded / (1024 * 1024)
                mb_tot = total_size / (1024 * 1024)
                sys.stdout.write(f"\r  Progress: {mb_down:.1f}/{mb_tot:.1f} MB ({percent:.1f}%)")
                sys.stdout.flush()
                
        elapsed = time.time() - start_time
        print(f"\n  [+] Transferred: {downloaded / (1024*1024):.2f} MB in {elapsed:.1f}s")


def download_and_extract():
    os.makedirs(DATA_DIR, exist_ok=True)
    
    for key, url in URLS.items():
        fname = os.path.basename(url)
        dest_path = os.path.join(DATA_DIR, fname)
        
        print(f"\n[*] Processing: {key} ({fname})")
        if not os.path.exists(dest_path) or os.path.getsize(dest_path) == 0:
            print(f"  Fetching from: {url}")
            download_file(url, dest_path)
        else:
            size_mb = os.path.getsize(dest_path) / (1024 * 1024)
            print(f"  [i] File already exists: {size_mb:.2f} MB")
            
        print(f"  Unpacking {fname} into {DATA_DIR}...")
        with tarfile.open(dest_path, "r:gz") as tar:
            tar.extractall(path=DATA_DIR)
        print(f"  [+] Extracted {key} successfully.")

    print("\n" + "="*60)
    print("VERIFICA INTEGRITA DATI ESTRATTI:")
    print("="*60)
    for root, dirs, files in os.walk(DATA_DIR):
        for f in files:
            if not f.endswith(".tar.gz"):
                fp = os.path.join(root, f)
                sz = os.path.getsize(fp)
                rel = os.path.relpath(fp, DATA_DIR)
                print(f"  - {rel:<45} : {sz:>12,} bytes")


if __name__ == "__main__":
    download_and_extract()
