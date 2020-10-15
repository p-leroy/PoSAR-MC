import json
import os.path
import re
import ipykernel
import requests
import logging

from urllib.parse import urljoin
from notebook.notebookapp import list_running_servers

def get_notebook_name():
    """
    Return the full path of the jupyter notebook.
    """
    kernel_id = re.search('kernel-(.*).json',
        ipykernel.connect.get_connection_file()).group(1)
    servers = list_running_servers()
    for ss in servers:
        response = requests.get(urljoin(ss['url'], 'api/sessions'),
            params={'token': ss.get('token', '')})
        for nn in json.loads(response.text):
            if nn['kernel']['id'] == kernel_id:
                relative_path = nn['notebook']['path']
                return os.path.join(ss['notebook_dir'], relative_path)

def loggingSetFh(name, logger, level=logging.INFO):
    # create file handler which logs even debug messages
    fh = logging.FileHandler(name, mode='w')
    fh.setLevel(level)
    # set the formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

def getNotebookLogger(level=logging.INFO):
    modulename = get_notebook_name().split('/')[-1].split('.ipynb')[0]
    logger = logging.getLogger()
    logger.setLevel(level)
    logfile = modulename + ".log"
    loggingSetFh(logfile, logger, level)
    return logger
