import hashlib
import os
import requests
import tempfile
from unittest import TestCase

from pyproteome import discoverer

PD_URL_BASE = (
    "https://media.githubusercontent.com/media/"
    "white-lab/pycamverter-data/master/"
)

PD14_URL_BASE = 'Test%20PD1.4/'
PD22_URL_BASE = 'Test%20PD2.2/'

PD14_URLS = [
    ("CK-C1-pY.msf", "96eda5b0e47f615cf000d1c5d3ecc8cd"),
    ("FAD-H1-Global.msf", "497ba1841faac883619d4f91a86a95cc"),
    ("FAD-H1-MPM2.msf", "5def08356bcfa5b679835e4a23dd1396"),
    ("Tau-4moHL1-Global.msf", "95d8089b5e4657b348bea0868c655478"),
]

PD22_URLS = [
    ('CK-H1-MPM2.msf', '0f880eb80b298233551a759864f25e44'),
]


class LoadSearchTest(TestCase):
    def fetch_url(self, url, md5hash):
        response = requests.get(url, stream=True)
        fd, path = tempfile.mkstemp(suffix='.msf')
        hash_md5 = hashlib.md5()

        with open(path, 'wb') as f:
            for block in response.iter_content(1024):
                hash_md5.update(block)
                f.write(block)

        if md5hash is not None and hash_md5.hexdigest() != md5hash:
            os.close(fd)
            raise Exception(
                "MD5 hash of {} does not match record: {} != {}"
                .format(url, md5hash, hash_md5.hexdigest())
            )

        return fd, path

    def test_load_pd14(self):
        for url, md5 in PD14_URLS:
            fd, path = self.fetch_url(
                PD_URL_BASE + PD14_URL_BASE + url,
                md5,
            )
            try:
                discoverer.read_discoverer_msf(url, msf_path=path)
            finally:
                os.close(fd)

    def test_load_pd22(self):
        for url, md5 in PD22_URLS:
            fd, path = self.fetch_url(
                PD_URL_BASE + PD22_URL_BASE + url,
                md5,
            )
            try:
                discoverer.read_discoverer_msf(url, msf_path=path)
            finally:
                os.close(fd)
