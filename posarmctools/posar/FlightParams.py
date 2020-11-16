import configparser, os

import logging
logger = logging.getLogger(__name__)

class FlightParams(object):
    def __init__(self, filename):
        self.filename = filename
        self.refs = self.__getReferencePoints()
        self.dir_sbg, self.session, self.sbg_hours = self.__getSbgParams()
        self.dir_posar, self.day, self.hours = self.__getRecordParams()
        
        self.root_dir = os.path.join(self.dir_posar, self.day)
        self.out_dir = os.path.join(self.dir_posar, self.day, "OUT")
        try:
            os.makedirs(self.out_dir, exist_ok=True) 
            logger.info(f"{self.out_dir} created successfully")
        except OSError as error:
            logger.error(f"{self.out_dir} can not be created, error: {error}")

    def __getReferencePoints(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [scene]
            refs = config.get("scene", "reference_points")
        return refs

    def __getSbgParams(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [sbg]
            dir_sbg = config.get("sbg", "dir_sbg")
            session = config.get("sbg", "session")
            sbg_hours = list(filter(None, config.get("sbg", "hours").splitlines()))
        return dir_sbg, session, sbg_hours

    def __getRecordParams(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [posar]
            dir_posar = config.get("posar", "dir_posar")
            day = config.get("posar", "day")
            lines = list(filter(None, config.get("posar", "hours").splitlines()))
            logger.info(f"dir_posar: {dir_posar}")
            logger.info(f"day: {day}")
            hours = {}
            for line in lines:
                if len(line.split(" ")) == 1:
                    hour = line
                    start = None
                    stop = None
                elif len(line.split(" ")) == 3:
                    hour = line.split(" ")[0]
                    start = line.split(" ")[1]
                    stop = line.split(" ")[2]
                    logger.info(f"trajectory information detected, start {start}, stop {stop}")
                hours[hour] = (start, stop)
                dirname = os.path.join(dir_posar, day, day + "_" + hour)
                if not os.path.isdir(dirname):
                    logger.error(f"{day}_{hour} does not exist")
                    break
                else:
                    logger.info(f"{day}_{hour} OK")
        return dir_posar, day, hours
