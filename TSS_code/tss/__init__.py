from tss.config import HOMER_PATH
import os

#/home/isshamie/software/homer/bin/
os.environ['PATH'] = HOMER_PATH + ':' + os.environ['PATH']
