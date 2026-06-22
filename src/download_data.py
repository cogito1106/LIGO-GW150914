import requests
import os

# LIGO data URLS from GWOSC for GW150914
# Two detectors: Hanford (H1) and Livingston (L1)
data_urls = {"H1": "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5",
    "L1": "https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3/L-L1_GWOSC_4KHZ_R1-1126259447-32.hdf5"
}

def download_file(url, destination):
    """Download a file from a URL to a local destination."""
    print(f"Downloading from {url} to {destination}...")
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Check if the request was successful

    with open(destination, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

if __name__ == "__main__":
    os.makedirs("data", exist_ok=True)  # Create a directory to store the data

    for detector, url in data_urls.items():
        filename = f"data/{detector}_GW150914.hdf5"
        download_file(url, filename)

    print("Done. Both detector data files have been downloaded.")
    