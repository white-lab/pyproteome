import os
import requests

from pyproteome import utils


def fetch_data(dirname, datas, base_url):
    utils.makedirs(dirname)

    for _, filename in datas.items():
        out_path = os.path.join(dirname, filename)

        if os.path.exists(out_path):
            continue

        url = base_url + filename
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(out_path, mode="wb") as f:
            for block in response.iter_content(1024):
                f.write(block)
