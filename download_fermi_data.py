#%%
import numpy as np
import requests
import time
from tqdm.contrib.concurrent import process_map
import os
#%%
week_in = 9
week_out = 745
gll_psc_v = 30
num_workers=None

download_photon_data=True
download_spacecraft_data=True
download_4FGL_catalogue=True

#%%
CBLUE = '\033[34m'
CEND = '\033[0m'

def print_msg_box(msg, indent=1, width=None, title=None):
    """Print message-box with optional title."""
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    print(CBLUE + box + CEND)


def download_url(args):
    t0 = time.time()
    url, fn = args[0], args[1]
    if os.path.exists(fn):
        # the file already exists, skipping
        return (url, time.time() - t0)
    else:
        try:
            r = requests.get(url)
            with open(fn, 'wb') as f:
                f.write(r.content)
            # print('url:', url, 'time (s):', time.time() - t0)
            return (url, time.time() - t0)
        except Exception as e:
            print('Exception in download_url():', e)



def download_parallel(iter, total, num_workers=None):
    if num_workers is None:
        num_workers = os.cpu_count()-1
    if total > num_workers*100:
        chuncksize = np.min([1, int(np.floor(total/(num_workers*100)))]) 
    else:
        chuncksize=1

    res = process_map(download_url, iter,
                      max_workers=num_workers, chunksize=chuncksize, total=total)
    return res

#%% download photon data
download_folder = "/data/users/Aurelio/fermi/photon"
os.makedirs(download_folder, exist_ok=True)

urls = []
destinations = []
for i in np.arange(week_in , week_out+1):
    urls.append(f"https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/photon/lat_photon_weekly_w{i:03d}_p305_v001.fits")
    destinations.append(f"{download_folder}/lat_photon_weekly_w{i:03d}_p305_v001.fits")

inputs = zip(urls, destinations)

#%%
if download_photon_data:
    print_msg_box("Downloading photon data")
    download_parallel(inputs, len(urls), num_workers=num_workers)
    print("done")
#%% download spacecraft
download_folder = "/data/users/Aurelio/fermi/spacecraft"
os.makedirs(download_folder, exist_ok=True)

url = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/lat_spacecraft_merged.fits"
destination = f"{download_folder}/lat_spacecraft_merged.fits"

inputs = (url, destination)
#%%
if download_spacecraft_data:
    print_msg_box("Downloading spacecraft file")
    download_url(inputs)
    print("done")
#%% download 4FGL sources
download_folder = "/data/users/Aurelio/fermi/4FGL"
os.makedirs(download_folder, exist_ok=True)
url = f"https://fermi.gsfc.nasa.gov/ssc/data/access/lat/12yr_catalog/gll_psc_v{gll_psc_v}.fit"
destination = f"{download_folder}/gll_psc_v{gll_psc_v}.fit"
inputs = (url, destination)
#%%
if download_4FGL_catalogue:
    print_msg_box("Downloading 4FGL catalogue")
    download_url(inputs)
    print("done")
#%%
print_msg_box("Done downloading data!")
