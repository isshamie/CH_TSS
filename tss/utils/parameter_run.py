import yaml
import glob
import os
import time
import pandas as pd
import numpy as np
from os.path import join


def wrap_mdir_params(dirs):
    if not (type(dirs) == list or type(dirs) == tuple):
        print("Needs to be a list or tuple")
        raise TypeError()
    created = list(map(lambda x: mdir_params(x), dirs))
    link_dirs(created)
    print(f"These directories were created: {created}")
    return created


def link_dirs(dirs):
    curr_str = "\n".join(dirs)
    for d in dirs:
        with open(join(d, ".linked_directories.txt"),"w") as f:
            f.write(curr_str)


def mdir_params(d):
    """
    Creates a directory if does not exist. If it does, creates the
    next numbered folder.
    e.g. mdir_params(foo/) -> foo/1 if no other numbered folder was
    there. if foo/ and foo/1 existed, foo/2 would be created. Time
    stamps and a hash function are created.
    :param d:
    :return:
    """
    if os.path.exists(d):
        # Get the numbered values
        count = 1
        new = os.path.join(d, str(count))
        while os.path.exists(new):
            count += 1
            new = join(d, str(count))
    else:
        new = join(d, "1")
    os.makedirs(new)

    time_file(new, name="DATE_created.txt")
    print(f"Creating directory {new}")

    return new


def mdir_old(d, replace=True, save_old=True, is_old=True):
    """ Function to make a directory. If the directory exists,
    it will either do nothing or replace the directory and save the
    old in a subfolder depending on the parameters
    Parameters:
        @d: type: str
            directory name.
        @replace: type=bool
                  default=True
                  To replace the old file or not
        @save_old: type=bool
                   default=True
                   Only on if replace=True. If true, will save the old
                   directory inside
                   'd/old_{index}' where the index will be the old
                   numbered version.

    Returns:
        True if the folder was not already present OR was
        replaced.
    """
    if os.path.exists(d):
        if replace:
            if save_old:
                if not is_old:
                    olds = glob.glob(os.path.join(d, "*"))
                    curr_index = len(olds)
                    count = 1
                    new_old = os.path.join(d, str(curr_index + count))
                    while os.path.exists(new_old):
                        count += 1
                        new_old = new_old[:new_old.rfind("_") + 1] + str(
                            count)
                else:
                    olds = glob.glob(os.path.join(d, "old_*"))
                    curr_index = len(olds)
                    count = 1
                    new_old = "old_%d" % (curr_index + count)
                    new_old = os.path.join(d, new_old)

                    # Keep incrementing until a new directory is found.
                    while os.path.exists(new_old):
                        count += 1
                        new_old = new_old[:new_old.rfind("_") + 1] + str(
                            count)

                os.mkdir(new_old)
                print("Directory already there. Moving contents into "
                      "%s" %
                      new_old)
                if is_old:
                    #This command moves everything but the old folder
                    cmd = "mv `ls -1 %s/* | grep -v old_` %s" % (d,
                                                                  new_old)
                    print(cmd)
                    os.system(cmd)

                time_file(new_old, name="DATE_moved.txt")

            else:
                cmd = "rm -rf %s/* " % d
                print("Removing old contents`")
                os.system(cmd)
        else:
            print("Directory already there. Keeping as is")
            return False
    else:
        os.makedirs(d)
        if not is_old:
            os.makedirs(os.path.join(d, "1"))
    return True


def time_file(outdir, name="DATE"):
    """ Function that generates a simple text file with a time stamp
    up to hours and minutes.
        @outdir: directory to save file in
        @name: file name to save to
        """
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y_%H%M', t)
    with open(os.path.join(outdir, name), 'w') as f:
        f.write(timestamp)

    # Random number to remember the same directory
    with open(os.path.join(outdir, name+"_randint"), 'w') as f:
        f.write("".join(np.random.randint(0, 10, 256).astype(str)))
    return


def append_name(f, prefix="", suffix="", to_rm=(".csv",".tsv")):
    """ Appends a prefix and suffix to string to the beginning of a
    filename
    that may be a full path """

    for i in to_rm:
        f = f.replace(i, "")

    if not prefix == "":
        prefix = prefix + "_"
    if not suffix == "":
        suffix = "_" + suffix

    f = os.path.join(os.path.dirname(f), prefix + os.path.basename(f)
                     + suffix)
    return f


def wrapper_yaml(yaml_dir):
    for param_f in glob.glob(yaml_dir):
        with open(param_f, 'r') as f:
            doc = yaml.load(f)
        #function = doc["function"]
        #cmd = "python {fnc} {f}"
    return


def run_function(d, func, create_dir=True):
    # for arg in sys.argv[1:]:
    #     print arg

    files = glob.glob(os.path.join(d, "*.yaml"))
    print(files)

    for param_f in files:
        with open(param_f, 'r') as f:
            params = yaml.load(f)
            new_d = os.path.join(params["prefix"], params["params"])

            if create_dir:
                # Make the directory
                mdir_params(new_d)
            run_function(new_d, func, params)
    return

