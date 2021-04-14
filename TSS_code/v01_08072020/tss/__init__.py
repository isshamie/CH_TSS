from tss.config import HOMER_PATH
from tss.data import *
import os

#/home/isshamie/software/homer/bin/
os.environ['PATH'] = HOMER_PATH + ':' + os.environ['PATH']
