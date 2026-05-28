from src import datalocation as dl
from src import runfile as rs

path = dl.fpath(r'.\data')
save = 'T'
show = 'F'
csv_save = 'T'


rs.runfile(path,save,show,csv_save)
