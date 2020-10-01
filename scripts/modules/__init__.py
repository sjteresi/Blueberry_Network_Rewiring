"""
Sundry helper functions.
"""


import logging
import errno
from functools import partial
from os import sysconf, strerror
import os

MAX_SYSTEM_RAM_GB = sysconf("SC_PAGE_SIZE") * sysconf("SC_PHYS_PAGES") / (1024.0 ** 3)
FILE_DNE = partial(FileNotFoundError, errno.ENOENT, strerror(errno.ENOENT))


def raise_if_no_file(filepath, logger=None, msg_fmt=None):

    logger = logger or logging.getLogger(__name__)
    msg_fmt = msg_fmt or "not a file:  %s"
    if not os.path.isfile(filepath):
        logger.critical(msg_fmt % filepath)
        raise FILE_DNE(filepath)


def raise_if_no_dir(filepath, logger=None, msg_fmt=None):

    logger = logger or logging.getLogger(__name__)
    msg_fmt = "not a directory:  %s"
    if not os.path.isdir(filepath):
        logger.critical(msg_fmt % filepath)
        raise FILE_DNE(filepath)
