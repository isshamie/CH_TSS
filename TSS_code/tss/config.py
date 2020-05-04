# -*- coding: utf-8 -*-
import click
import logging
import yaml
from dotenv.main import dotenv_values
import pathlib
import os

#HOMER_PATH = "/home/isshamie/software/homer/bin/"
HOMER_PATH = "/data2/resources/software/homer/v4.11_10242019/bin//"

## Get the project directory
def set_variables(params):
    # Paths
    global PREFIX
    PREFIX = params["PREFIX"]
    global ROOT_DIR
    ROOT_DIR = params["ROOT_DIR"]
    global DATA_DIR
    DATA_DIR = params["DATA_DIR"]
    global DATA_PROCESSED_DIR
    DATA_PROCESSED_DIR = params["DATA_PROCESSED_DIR"]
    global SUPPLEMENTAL_DIR
    SUPPLEMENTAL_DIR = params["SUPPLEMENTAL_DIR"]
    global PIPELINE
    PIPELINE = params["PIPELINE"]
    ## Genome Files
    global GENOME_DIR
    GENOME_DIR = params["GENOME_DIR"]
    global GENOME_FA
    GENOME_FA = params["GENOME_FA"]
    global GENOME_GFF3
    GENOME_GFF3 = params["GENOME_GFF3"]
    global GENOME_GTF
    GENOME_GTF = params["GENOME_GTF"]
    ## Params
    global TISSUES
    TISSUES = params["TISSUES"]
    # meta file for names
    global META_FILE
    META_FILE = params["META_FILE"]
    return



@click.command()
@click.argument('parameter_file', type=click.Path(exists=True))
def main(parameter_file):
    """ Runs data processing pipeline to turn tags and peaks data
    from DATA_DIR to DATA_PROCESSED_DIR.
    """
    logger = logging.getLogger(__name__)
    logger.info('making final data set from raw data')
    PARAMS = yaml.load(open(parameter_file), Loader=yaml.FullLoader)
    run(PARAMS)


def run(params):
    set_variables(params)
    return

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    main()

else:
    PARAMS = dotenv_values()
    print("params", PARAMS)
    #run(PARAMS)
