import logging

import pytest
import sys
import os
import zipfile


import requests
from tqdm import tqdm


def download_from_url(url, dst, desc=None, overwrite=False):
    file_size = int(requests.head(url).headers["Content-Length"])
    if overwrite:
        logging.info(f"overwriting {dst}")
        first_byte = 0
    elif os.path.exists(dst):
        first_byte = os.path.getsize(dst)
    else:
        first_byte = 0
    if first_byte >= file_size:
        return file_size
    header = {"Range": "bytes=%s-%s" % (first_byte, file_size)}
    if desc is None:
        desc = url.split("/")[-1]
    pbar = tqdm(
        total=file_size, initial=first_byte, unit="B", unit_scale=True, desc=desc
    )
    req = requests.get(url, headers=header, stream=True)
    with (open(dst, "ab")) as f:
        for chunk in req.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
                pbar.update(1024)
    pbar.close()

    return file_size


def pytest_sessionstart(session):
    """
    Called after the Session object has been created and
    before performing collection and entering the run test loop.
    """
    zip_dest = os.path.join(os.path.dirname(__file__), "data.zip")
    dir_dest = os.path.join(os.path.dirname(__file__), "data")

    if os.path.isdir(dir_dest):
        logging.info(
            f"Test data found. Not downloading. To force download, delete {dir_dest}"
        )
    else:
        logging.info("no test data found. downloading.")

        download_from_url(
            "https://ucc28400d85c5cc3ad392656f954.dl.dropboxusercontent.com/zip_download_get/Aedl8JwZHiq_trPNj0ad27K2Q4qQFqOJUalz6RmDS_TpbOYerd4giEyHQSnAL7pEHPDP2p8UqLXT_v89XWfoYkgiHtb47BflKfL9kUvr74pW6Q",
            zip_dest,
            overwrite=False,
            desc="Downloading Test Data",
        )
        with zipfile.ZipFile(zip_dest, "r") as zip_ref:
            zip_ref.extractall(dir_dest)
        os.remove(zip_dest)


def pytest_load_initial_conftests(args):
    if "xdist" in sys.modules:  # pytest-xdist plugin
        import multiprocessing

        num = max(multiprocessing.cpu_count() / 2, 1)
        args[:] = ["-n", str(num)] + args


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
