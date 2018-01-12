import os
import requests


def fetch_data(dir, datas, base_url):
    try:
        os.makedirs(dir)
    except:
        pass

    for _, filename in datas.items():
        out_path = os.path.join(dir, filename)

        if os.path.exists(out_path):
            continue

        url = base_url + filename
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(out_path, mode="wb") as f:
            for block in response.iter_content(1024):
                f.write(block)
