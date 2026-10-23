#!/usr/bin/env python3
"""
download_cohort.py
==================
Download and verify the multi-cancer cohort (Breast Cancer & Melanoma)
from authentic 10x Genomics CytAssist FFPE repositories.
"""

import os
import sys
import tarfile
import urllib.request
import time
import subprocess

BASE_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

DATASETS = {
    "breast": {
        "spatial": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz",
        "matrix": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.tar.gz"
    },
    "melanoma": {
        "spatial": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz",
        "matrix": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Skin_Melanoma/CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.tar.gz"
    }
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
                sys.stdout.write(f"\r  Downloading: {mb_down:.1f}/{mb_tot:.1f} MB ({percent:.1f}%)")
                sys.stdout.flush()
                
        elapsed = time.time() - start_time
        rate = (downloaded / (1024 * 1024)) / max(elapsed, 0.01)
        print(f"\n  [+] Downloaded: {downloaded / (1024*1024):.2f} MB in {elapsed:.1f}s ({rate:.2f} MB/s)")


def verify_gzip_integrity(file_path):
    res = subprocess.run(["gzip", "-t", file_path], capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(f"GZIP CRC/integrity error on {file_path}: {res.stderr.decode()}")
    print(f"  [✓] Gzip integrity verified (RFC 1952 compliant): {os.path.basename(file_path)}")


def process_dataset(name, urls):
    target_dir = os.path.join(BASE_DATA_DIR, name)
    os.makedirs(target_dir, exist_ok=True)
    print(f"\n{'='*60}\nINIZIO INGESTION: {name.upper()}\n{'='*60}")
    
    for key, url in urls.items():
        fname = os.path.basename(url)
        tar_path = os.path.join(target_dir, fname)
        
        print(f"\n[*] Componente: {key} ({fname})")
        if not os.path.exists(tar_path) or os.path.getsize(tar_path) == 0:
            print(f"  Fetching: {url}")
            download_file(url, tar_path)
        else:
            print(f"  [i] Archivio già presente: {os.path.getsize(tar_path)/(1024*1024):.2f} MB")
            
        verify_gzip_integrity(tar_path)
        
        print(f"  Estrazione in {target_dir}...")
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=target_dir)
        print(f"  [+] Estrazione {key} completata.")


def audit_cohort():
    print(f"\n{'='*60}\nAUDIT FINALE DELLA COORTE MULTI-CANCRO (10x GENOMICS)\n{'='*60}")
    for tumor in ["lung", "breast", "melanoma"]:
        tdir = os.path.join(BASE_DATA_DIR, tumor)
        print(f"\n--- Tumore: {tumor.upper()} ({tdir}) ---")
        if not os.path.exists(tdir):
            print("  [!] Cartella mancante")
            continue
        for root, dirs, files in os.walk(tdir):
            for f in sorted(files):
                if not f.endswith(".tar.gz"):
                    fp = os.path.join(root, f)
                    sz = os.path.getsize(fp)
                    rel = os.path.relpath(fp, tdir)
                    print(f"  - {rel:<45} : {sz:>12,} bytes")


if __name__ == "__main__":
    for tumor_name, url_dict in DATASETS.items():
        process_dataset(tumor_name, url_dict)
    audit_cohort()
