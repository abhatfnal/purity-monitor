from os import listdir
from os.path import isfile, join
import numpy as np
from optparse import OptionParser

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--filepath", action="store", type="string", dest="filepath", default="", help="Specify the path.")
(options, args) = parser.parse_args()

onlyfiles = [f for f in listdir(options.filepath) if isfile(join(options.filepath, f))]

print onlyfiles
