import sys,argparse,csv

parser = argparse.ArgumentParser(
    prog='python3 evolution_analyzer.py',
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input-files",help="Input files to be investigated. Program accepts any number of input files.",nargs='+')
parser.add_argument("-r","--reference",help="Reference file for sequences investigation.")
parser.add_argument("-t","--top",type=int,help="Number of investigated sequences. Default is 10.",default=10)
parser.add_argument("-o","--output",help="Output File. Mandatory.")
args = parser.parse_args()

if args.output == '':
    print("No output file has been defined. Please define an output file before launching Evolution Analyzer.")
    sys.exit()

Ref = {}
header = []
header.append(str(args.reference))
with open(args.reference,'r') as file:    #create the reference lists for investigated sequences
 reference = file.readlines()          #and relative frequencies (ref_seq and ref_freq)
 for i in range(args.top):
  Ref[reference[i].rstrip('\n').split(':')[0]] = [float(reference[i].rstrip('\n').split(':')[1])]

for item in args.input_files:             #read every input file and temporary stores
 header.append(str(item))
 temp = {}
 with open(item,'r') as file:          #sequences and frequencies
  for line in file:
   temp[line.strip('\n').split(':')[0]] = float(line.strip('\n').split(':')[1]) 
  comm = list(set(Ref.keys()).intersection(set(temp.keys())))
  for seq in Ref.keys():
   if seq in comm:
    Ref[seq].append(temp[seq])
   else:
    Ref[seq].append('Not Found')
  
w = csv.writer(open(args.output,'w'), dialect='excel')
w.writerow(header)
for obj in Ref.keys():
 w.writerow([obj+','+','.join(str(x) for x in Ref[obj])])
