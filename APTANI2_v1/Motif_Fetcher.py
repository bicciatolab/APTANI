import sys,argparse,csv
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
    prog='python3 Motif_Fetcher.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=("""Motif Fetcher processes an APTANI² output and retrieves the motifs and score values associated to the investigated aptamer sequences."""))
parser.add_argument("-v","--version", action="version",version="Motif Fetcher V1.0 developed by Jimmy Caroli\nBicciato's Bioinformatics Laboratory\nUniversity of Modena and Reggio Emilia 2017")
parser.add_argument("-a","--aptamers",help="Aptamers to be investigated. Program accepts any number of sequences files.",nargs='+')
parser.add_argument("-i","--filein",type=str, help="Input file if different from Motifs.csv", default ="")
args = parser.parse_args()

if args.filein == "":
 Input = pd.read_csv("Motifs.csv",sep=" ")
else:
 Input = pd.read_csv(args.filein,sep=" ")

for item in args.aptamers:
  w = csv.writer(open('Motifs_in_{0}.csv'.format(item), 'w'))
  Indexes = Input.index[Input['Aptamer'] == item].tolist()
  for number in Indexes:
    w.writerow([Input.iloc[number]['Motif'] +': '+ str(Input.iloc[number]['Score'])])

