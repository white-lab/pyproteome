import os
import requests

import pyproteome as pyp


def fetch_data(dirname, datas, base_url):
    pyp.utils.makedirs(dirname)

    for name in datas:
        filename = name + ".msf"
        out_path = os.path.join(dirname, filename)

        if os.path.exists(out_path):
            continue

        url = base_url + filename
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(out_path, mode="wb") as f:
            for block in response.iter_content(1024):
                f.write(block)
