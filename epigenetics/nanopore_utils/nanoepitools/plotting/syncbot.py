from pathlib import Path
import yaml
import subprocess
from typing import Union

import logging


class SyncBot:
    """
    Helps sync plots to dash. Assumes that "dash" is configured in SSH config.
    """
    
    def __init__(self, configfile: Union[str, Path] = None):
        self.logger = logging.getLogger("syncbot")
        handler = logging.StreamHandler()
        formatter = logging.Formatter("%(asctime)s [%(name)s, %(levelname)s] - %(message)s")
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        
        if configfile is None:
            configfile = Path.home().joinpath(".config/syncplots.yaml")
        
        if not isinstance(configfile, Path):
            configfile = Path(configfile)
        
        if configfile.is_file():
            with open(configfile, "r") as f:
                self.config = yaml.safe_load(f)
        else:
            self.logger.warning("No nanoepitools syncbot config. Will not synchronize plots.")
            self.config = {}
    
    def sync(self):
        for key in self.config:
            subconf = self.config[key]
            rsync = f"rsync -avhL --delete {subconf['path']} dash:/data/syncplots/{subconf['name_on_server']}"
            if "exclude" in subconf:
                rsync += " --exclude=".join([""] + subconf["exclude"])
            
            self.logger.debug(f"syncbot running command: {rsync}")
            process = subprocess.Popen(rsync.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with process.stdout:
                for line in process.stdout.readlines():
                    self.logger.debug(line.decode().strip())
            process.wait()
