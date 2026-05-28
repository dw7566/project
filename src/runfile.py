from .tocsv import csv_mod
from .ivfitting import IV

import glob


def runfile(path, save, show, csv_save):
    xml = []
    for filename in glob.glob(path, recursive=True):
        xml.append(filename)

    print(xml)
    for i in xml:
        IV(i,save,show)
        csv_mod(i,csv_save)
