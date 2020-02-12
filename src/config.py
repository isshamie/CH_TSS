# -*- coding: utf-8 -*-
import click
import logging
import yaml

import pathlib
import os







## Get the project directory
def run():
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

    # Paths
    PREFIX = PARAMS["PREFIX"]
    ROOT_DIR = PARAMS["ROOT_DIR"]
    DATA_DIR = PARAMS["DATA_DIR"]
    DATA_PROCESSED_DIR = PARAMS["DATA_PROCESSED_DIR"]
    SUPPLEMENTAL_DIR = PARAMS["SUPPLEMENTAL_DIR"]

    PIPELINE = PARAMS["PIPELINE"]
    ## Genome Files
    GENOME_DIR = PARAMS["GENOME_DIR"]
    GENOME_FA = PARAMS["GENOME_FA"]
    GENOME_GFF3 = PARAMS["GENOME_GFF3"]
    GENOME_GTF = PARAMS["GENOME_GTF"]
    ## Params
    TISSUES = PARAMS["META_FILE"]
    # meta file for names
    META_FILE = PARAMS["META_FILE"]
    run()



if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    main()

else:
    ## Get the project directory
    path = os.path.abspath(__file__)
    dir_path = pathlib.Path(os.path.dirname(path))
    ROOT_DIR = dir_path.parents[0]
    print(ROOT_DIR)


    parameter_file = "parameters/params.yaml"
    PARAMS = yaml.load(open(parameter_file), Loader=yaml.FullLoader)
    # Paths
    PREFIX = PARAMS["PREFIX"]
    ROOT_DIR = PARAMS["ROOT_DIR"]
    DATA_DIR = PARAMS["DATA_DIR"]
    DATA_PROCESSED_DIR = PARAMS["DATA_PROCESSED_DIR"]
    SUPPLEMENTAL_DIR = PARAMS["SUPPLEMENTAL_DIR"]

    PIPELINE = PARAMS["PIPELINE"]
    ## Genome Files
    GENOME_DIR = PARAMS["GENOME_DIR"]
    GENOME_FA = PARAMS["GENOME_FA"]
    GENOME_GFF3 = PARAMS["GENOME_GFF3"]
    GENOME_GTF = PARAMS["GENOME_GTF"]
    ## Params
    TISSUES = PARAMS["META_FILE"]
    # meta file for names
    META_FILE = PARAMS["META_FILE"]

