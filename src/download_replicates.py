#!/usr/bin/env python3
"""
download_replicates.py
======================
Download authentic 10x Genomics Visium v1.1.0 consecutive serial sections
(V1_Breast_Cancer_Block_A_Section_1 and Section_2) to evaluate technical and biological reproducibility.
"""

import os
import sys
import tarfile
import urllib.request
import time
import subprocess

BASE_DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

DATASETS = {
    "breast_visium_sec1": {
        "spatial": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz",
        "matrix": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.tar.gz"
    },
    "breast_visium_sec2": {
        "spatial": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_spatial.tar.gz",
        "matrix": "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.tar.gz"
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
        chunk_size = 1024 * 1024
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
        print(f"\n  [+] Transferred: {downloaded / (1024*1024):.2f} MB in {elapsed:.1f}s")


def verify_gzip(path):
    res = subprocess.run(["gzip", "-t", path], capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(f"GZIP corrupted: {path}")
    print(f"  [✓] RFC 1952 Gzip Verified: {os.path.basename(path)}")


def process_all():
    for name, urls in DATASETS.items():
        target_dir = os.path.join(BASE_DATA_DIR, name)
        os.makedirs(target_dir, exist_ok=True)
        print(f"\n{'='*60}\nINGESTING SERIAL REPLICATE: {name.upper()}\n{'='*60}")
        for k, u in urls.items():
            fn = os.path.basename(u)
            dest = os.path.join(target_dir, fn)
            if not os.path.exists(dest) or os.path.getsize(dest) == 0:
                print(f"[*] Downloading {k} from {u}...")
                download_file(u, dest)
            verify_gzip(dest)
            print(f"[*] Extracting {fn} to {target_dir}...")
            with tarfile.open(dest, "r:gz") as tar:
                tar.extractall(path=target_dir)
            print(f"  [+] Extracted {k}.")


if __name__ == "__main__":
    process_all()
