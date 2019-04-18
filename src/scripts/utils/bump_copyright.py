#!/usr/bin/env python
"""
Walk the directory tree and replace the '20XX' copyright dates in lines with
'Copyright 2003-20XX' with '2016'. Usage:

    ./bump_copyright.py          # use "." as the start directory
    ./bump_copyright.py dirname  # use dirname as the start directory
"""

from __future__ import print_function
import os
import re
import sys


search_string = r'Copyright \(c\) 2003-20\d\d'
replace_string = r'Copyright (c) 2003-2019'


def has_hidden_dirs(dirpath):
    """
    Check if the directory path has any directories leading with a '.' (accept
    the current working directory as okay).
    """
    dirs = dirpath.split("/")
    for dir in dirs:
        # match - anchor a "." at the beginning, then 1 or more of any char
        m = re.match(r'^\.(.+)', dir)
        if m:
            return True
    return False


def get_the_paths(root_dir):
    """
    Get all the files NOT in hidden directories
    """
    path_collection = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        use_dir = not has_hidden_dirs(dirpath)
        if use_dir:
            for filename in filenames:
                # Don't try to replace the strings in this script
                if not re.search(sys.argv[0].split("/")[-1], filename):
                    fullpath = os.path.join(dirpath, filename)
                    path_collection.append(fullpath)
    return path_collection


def do_the_subs(path_collection, search, replace):
    """
    Look for the copyright string and sub the date
    """
    for path in path_collection:
        with open(path, 'r+') as p:
            contents = p.read()
            pattern = re.compile(search)
            contents = pattern.sub(replace, contents)
            p.seek(0)
            p.truncate()
            p.write(contents)


if __name__ == '__main__':

    if '-h' in sys.argv or '--help' in sys.argv:
        print(__doc__)
        sys.exit(1)

    d = "."
    if len(sys.argv) > 1:
        d = sys.argv[1]

    path_collection = get_the_paths(d)
    do_the_subs(path_collection, search_string, replace_string)
