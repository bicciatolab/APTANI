def percentage(p, w):
  return(int((p*w)/100))

def word_count(filename,num,frequency):
  words = []
  subopt={}
#  pattern = re.compile('[ATGC]{%s}'%num)
  with open(filename,'r') as f:
    for i in islice(f,1,None,4):
      if len(i.strip()) == num:
        words.append(i.strip())
#      if pattern.match(i):
#        words.append(i.strip())
  tot = float(len(words))
  counter = collections.Counter(words)
  counter = dict(counter)
  counter = {k: v / float(len(words)) for k, v in counter.items()}

  lista=sorted(counter.items(),key=lambda x:x[1],reverse=True)

  w = csv.writer(open('count.txt', 'w'))
  for (word,freq) in lista:
    w.writerow([str(word)+': '+str(freq)])
    if freq > frequency:
      subopt[word]=freq
  
  w = csv.writer(open('subopt_input.txt', 'w'))
  for key, val in subopt.items():
    w.writerow([key])

def word_major(filename,num,frequency):
  words = []
  subopt={}
  with open(filename,'r') as f:
    for i in islice(f,1,None,4):
      if len(i.strip()) >= num:
        words.append(i.strip())
  
  tot = float(len(words))
  counter = collections.Counter(words)
  counter = dict(counter)
  counter = {k: v / float(len(words)) for k, v in counter.items()}

  lista=sorted(counter.items(),key=lambda x:x[1],reverse=True)

  w = csv.writer(open('count.txt', 'w'))
  for (word,freq) in lista:
    w.writerow([str(word)+': '+str(freq)])
    if freq > frequency:
      subopt[word]=freq

  w = csv.writer(open('subopt_input.txt', 'w'))
  for key, val in subopt.items():
    w.writerow([key])

def variable_region(filename,frequency,cutoff):
  words = []
  subopt = {}
  count = {}
  with open(filename,'r') as f:
    for i in islice(f,1,None,4):
      try:
        len(re.search(r'{0}(.*?){1}'.format(args.left.replace('U','T'),args.right.replace('U','T')), i).group(1)) > cutoff and words.append(re.search(r'{0}(.*?){1}'.format(args.left.replace('U','T'),args.right.replace('U','T')), i).group(1))
      except AttributeError:
        continue

  tot = float(len(words))
  counter = collections.Counter(words)
  counter = dict(counter)
  counter = {k: v / float(len(words)) for k, v in counter.items()}

  lista=sorted(counter.items(),key=lambda x:x[1],reverse=True)

  w = csv.writer(open('count.txt', 'w'))
  for (word,freq) in lista:
    w.writerow([str(word)+': '+str(freq)])
    if freq > frequency:
      subopt[args.struct_left.replace('U','T')+word+args.struct_right.replace('U','T')]=freq

  w = csv.writer(open('subopt_input.txt', 'w'))
  for key, val in subopt.items():
    w.writerow([key])

def hairpin(DATA,MFE,Aptamers,Score,StructScore):
  hp = re.compile(r"\([^()]+\)")
  tot = 0                                 # counter per i motivi
  Dict = {}
  hp = re.compile(r"\([^()]+\)")
  tot = 0                                 # counter per i motivi
  Dict = {}
  for k in DATA.keys():             # per ogni aptamero
    if (LEFT == '') and (RIGHT == ''):
      LimitL = 0
      LimitR = len(k)
    else:
      LimitL = len(LEFT)
      LimitR = len(k)-len(RIGHT)
    for x in DATA[k]:               # per ogni struttura secondaria calcolata
      for m in hp.finditer(x[0]):          # trovo le hairpin nella struttura secondaria
        if (m.end() <= LimitL) or (m.start() >= LimitR):
          continue
        else: 
          tot += 1                           # conteggio le hairpin
          StructScore[k][x[0]].append([k[m.start():m.end()], str(m.start()+1)+'-'+str(m.end())])
          if k[m.start():m.end()] not in Aptamers[k]:
            Aptamers[k].append(str(k[m.start():m.end()]))
          if k[m.start():m.end()] in Dict:   # se ho già trovato l'hairpin
            Dict[k[m.start():m.end()]][0] += float(x[-1])    # sommo il valore di energia
            Dict[k[m.start():m.end()]][1] += 1               # conteggio quante volte trovo questa hairpin
          else:                              # se hairpin e' nuova
            Dict[k[m.start():m.end()]] = [float(x[-1]),1]    # salvo nuova chiave nel dizionario

  for key in Dict:
    Score[key] = (Dict[key][1]/tot)**(1-((Dict[key][0]/Dict[key][1])/MFE))

  header = 'Motif Score'
  w = csv.writer(open('Hairpin_Score.txt', 'w'))
  w.writerow([header])
  for key in Score:
    w.writerow([str(key)+' '+str(Score[key])])

def quadruplex(DATA,MFE,Aptamers,Score,Plex,StructScore):
  quad = re.compile(r'(\+{2,}.{1,7}){3}\+{2,}')
  tot = 0
  Dict = {}

  for k in DATA.keys():
    if (LEFT == '') and (RIGHT == ''):
      LimitL = 0
      LimitR = len(k)
    else:
      LimitL = len(LEFT)
      LimitR = len(k) - len(RIGHT)

    for x in DATA[k]:               # per ogni struttura secondaria calcolata
      for m in quad.finditer(x[0]):          # trovo le hairpin nella struttura secondaria
        if (m.end() <= LimitL) or (m.start() >= LimitR):
          continue
        else:
          tot += 1                           # conteggio le hairpin
          if k[m.start():m.end()] not in Aptamers[k]:
            Aptamers[k].append(str(k[m.start():m.end()]))
            Plex.append(k)
            StructScore[k][x[0]].append([k[m.start():m.end()], str(m.start()+1)+'-'+str(m.end())])
          if k[m.start():m.end()] in Dict:   # se ho già trovato l'hairpin
            Dict[k[m.start():m.end()]][0] += float(x[-1])    # sommo il valore di energia
            Dict[k[m.start():m.end()]][1] += 1               # conteggio quante volte trovo questa hairpin
            StructScore[k][x[0]].append([k[m.start():m.end()], str(m.start()+1)+'-'+str(m.end())])
          else:                              # se hairpin e' nuova
            Dict[k[m.start():m.end()]] = [float(x[-1]),1]    # salvo nuova chiave nel dizionario
            StructScore[k][x[0]].append([k[m.start():m.end()], str(m.start()+1)+'-'+str(m.end())])

  for key in Dict:
    Score[key] = (Dict[key][1]/tot)**(1-((Dict[key][0]/Dict[key][1])/MFE))

  header = 'Motif Score'
  w = csv.writer(open('Quadruplex_Score.txt', 'w'))
  w.writerow([header])
  for key in Score:
    w.writerow([str(key)+' '+str(Score[key])])

def right_bulges(DATA,MFE,Aptamers,Score,StructScore):
  left = re.compile(r"(\(\()")
  right = re.compile(r"(\)\.+\))")
  abort = re.compile(r"(\.+\))")
  bulges = {}
  tot = 0

  for key in DATA.keys():
    if (LEFT == '') and (RIGHT == ''):
      LimitL = 0
      LimitR = len(key)
    else:
      LimitL = len(LEFT)
      LimitR = len(key) - len(RIGHT)  
    for item in DATA[key]:
      for m in left.finditer(item[0]):
        for l in right.finditer(item[0][m.end():]):
          if (item[0][m.end():l.start()+m.end()].count('(') == item[0][m.end():l.start()+m.end()].count(')')) and (item[0][m.end()+1] != '.') and ((m.start() in range(LimitL,LimitR)) or (m.end() in range(LimitL,LimitR)) or (l.start()+m.end() in range(LimitL,LimitR)) or (l.end()+m.end() in range(LimitL,LimitR))):
            tot += 1 
            StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(m.start())+'-'+str(m.end()),str(l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
            if str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]) not in Aptamers[key]:
              Aptamers[key].append(str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]))
            if str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]) in bulges:
              bulges[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])][0] += float(item[-1])
              bulges[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])][1] += 1
              StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(1+m.start())+'-'+str(m.end()),str(1+l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
              break
          else:
            bulges[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])] = [float(item[-1]),1]
            StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(1+m.start())+'-'+str(m.end()),str(1+l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
            break

  for key in bulges:
    Score[key] = (bulges[key][1]/tot)**(1-((bulges[key][0]/bulges[key][1])/MFE))

  header = 'Motif Score'
  w = csv.writer(open('Right_Bulges_Score.txt', 'w'))
  w.writerow([header])
  for key in Score:
    w.writerow([str(key)+' '+str(Score[key])])

def left_bulges(DATA,MFE,Aptamers,Score,StructScore):
  left = re.compile(r"(\(\.+\()")
  right = re.compile(r"(\)\))")
  abort = re.compile(r"(\.+\))")
  bulges = {}
  tot = 0

  for key in DATA.keys():
    if (LEFT == '') and (RIGHT == ''):
      LimitL = 0
      LimitR = len(key)
    else:
      LimitL = len(LEFT)
      LimitR = len(key) - len(RIGHT)  
    for item in DATA[key]:
      for m in left.finditer(item[0]):
        for l in right.finditer(item[0][m.end():]):
          if (item[0][m.end():l.start()+m.end()].count('(') == item[0][m.end():l.start()+m.end()].count(')')) and (item[0][l.start()+m.end()-1] != '.') and ((m.start() in range(LimitL,LimitR)) or (m.end() in range(LimitL,LimitR)) or (l.start()+m.end() in range(LimitL,LimitR)) or (l.end()+m.end() in range(LimitL,LimitR))): 
            tot += 1
            StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(m.start())+'-'+str(m.end()),str(l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
            if str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]) not in Aptamers[key]:
              Aptamers[key].append(str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]))
            if str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]) in bulges:
              bulges[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])][0] += float(item[-1])
              bulges[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])][1] += 1
              StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(1+m.start())+'-'+str(m.end()),str(1+l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
              break
            else:
              bulges[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])] = [float(item[-1]),1]
              StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(1+m.start())+'-'+str(m.end()),str(1+l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
              break

  for key in bulges:
    Score[key] = (bulges[key][1]/tot)**(1-((bulges[key][0]/bulges[key][1])/MFE))

  header = 'Motif Score'
  w = csv.writer(open('Left_Bulges_Score.txt', 'w'))
  w.writerow([header])
  for key in Score:
    w.writerow([str(key)+' '+str(Score[key])])

def intra(DATA,MFE,Aptamers,Score,StructScore):
  left = re.compile(r"(\(\.+\()")
  right = re.compile(r"(\)\.+\))")
  abort = re.compile(r"(\.+\))")
  intra = {}
  tot = 0

  for key in DATA.keys():
    if (LEFT == '') and (RIGHT == ''):
      LimitL = 0
      LimitR = len(key)
    else:
      LimitL = len(LEFT)
      LimitR = len(key) - len(RIGHT)  
    for item in DATA[key]:
      for m in left.finditer(item[0]):
        for l in right.finditer(item[0][m.end():]):
          if (item[0][m.end():l.start()+m.end()].count('(') == item[0][m.end():l.start()+m.end()].count(')')) and ((m.start() in range(LimitL,LimitR)) or (m.end() in range(LimitL,LimitR)) or (l.start()+m.end() in range(LimitL,LimitR)) or (l.end()+m.end() in range(LimitL,LimitR))): 
            tot += 1
            StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(m.start())+'-'+str(m.end()),str(l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
            if str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]) not in Aptamers[key]:
              Aptamers[key].append(str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]))
            if str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()]) in intra:
              intra[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])][0] += float(item[-1])
              intra[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])][1] += 1
              StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(1+m.start())+'-'+str(m.end()),str(1+l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
              break
            else:
              intra[str(key[m.start():m.end()])+str('|')+str(key[m.end()+l.start():m.end()+l.end()])] = [float(item[-1]),1]
              StructScore[key][item[0]].append([str(key[m.start():m.end()])+'|'+str(key[m.end()+l.start():m.end()+l.end()]),str(1+m.start())+'-'+str(m.end()),str(1+l.start()+m.end())+'-'+str(l.end()+m.end()+1)])
              break

  for key in intra:
    Score[key] = (intra[key][1]/tot)**(1-((intra[key][0]/intra[key][1])/MFE))
 
  header = 'Motif Score'
  w = csv.writer(open('Intra_Score.txt', 'w'))
  w.writerow([header])
  for key in Score:
    w.writerow([str(key)+' '+str(Score[key])])

import sys,math,random,re,subprocess,string,locale,os,time,argparse,collections,csv
from itertools import islice
from collections import Counter

parser = argparse.ArgumentParser(
    prog='python3.4 APTANI2.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=("""

 .----------------.  .----------------.  .----------------.  .----------------.  .-----------------. .----------------.  .----------------. 
| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |
| |      __      | || |   ______     | || |  _________   | || |      __      | || | ____  _____  | || |     _____    | || |    _____     | |
| |     /  \     | || |  |_   __ \   | || | |  _   _  |  | || |     /  \     | || ||_   \|_   _| | || |    |_   _|   | || |   / ___ `.   | |
| |    / /\ \    | || |    | |__) |  | || | |_/ | | \_|  | || |    / /\ \    | || |  |   \ | |   | || |      | |     | || |  |_/___) |   | |
| |   / ____ \   | || |    |  ___/   | || |     | |      | || |   / ____ \   | || |  | |\ \| |   | || |      | |     | || |   .'____.'   | |
| | _/ /    \ \_ | || |   _| |_      | || |    _| |_     | || | _/ /    \ \_ | || | _| |_\   |_  | || |     _| |_    | || |  / /____     | |
| ||____|  |____|| || |  |_____|     | || |   |_____|    | || ||____|  |____|| || ||_____|\____| | || |    |_____|   | || |  |_______|   | |
| |              | || |              | || |              | || |              | || |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------' 

APTANI² processes a SELEX or Cell-SELEX output and retrieves the most probable aptamer-related binding motifs"""))
parser.add_argument("-v","--version", action="version",version="APTANI² V2.0 developed by Jimmy Caroli\nBicciato's Bioinformatics Laboratory\nUniversity of Modena and Reggio Emilia 2018")
parser.add_argument("Input File", help="SELEX/Cell-SELEX output data to be analyzed. Positional argument at the end of the command line. Format accepted is FastQ")
parser.add_argument("-e","--energy", type=int, help="Energy Range in RNAsubopt calculation. Default is 3", default=3)
parser.add_argument("-f","--frequency",type=float, help="Frequency cutoff for aptamers processing. Default is 0.000001", default=0.000001)
parser.add_argument("-n","--number",type=int,help="Lenght of the Aptamer sequence in SELEX/Cell-SELEX output. Default is 99", default=99)
parser.add_argument("-m","--major",action="store_true",help="Select if input file has sequences of variable lengths. Parameter -n will define the minimum length of the sequence to be searched")
parser.add_argument("-u","--upper",type=float,help="Upper percentage of the distribution of the MtfScore values used to define the top scoring motifs. Default if to select first percentile (Top 1%%) of scored motifs.", default=1) 
parser.add_argument("--left",type=str, help="Left Tag to be searched. Input as RNA nt [AUCG]. Default is AUGCGG", default= 'AUGCGG')
parser.add_argument("-L","--struct-left", type=str, help="Left tag to be attached during secondary structure calculations. Input as RNA nt [AUCG]. Type '-L NONE' if you do not want a left structural tag. No deafault is set.", default= '')
parser.add_argument("-R","--struct-right", type=str, help="Right tag to be attached during secondary structure calculations. Input as RNA nt [AUCG]. Type '-R NONE' if you do not want a right structural tag. No default is set.", default='')
parser.add_argument("--right",type=str, help="Right Tag to be searched. Input as RNA nt [AUCG]. Default is CAGACG", default='CAGACG')
parser.add_argument("-x","--variable",action="store_true",help="Searches and elaborates ONLY the variable region of aptamers comprised between the defined tags")
parser.add_argument("--cutoff",type=int,help="Lenght cutoff for the searched variable region between left and right tags. Default is 30",default=30)
args = parser.parse_args()

filename = sys.argv[-1]

NONE = ['NONE','None','none','no','NO','']

if ((args.struct_left in NONE) and (args.struct_right not in NONE)) or ((args.struct_left not in NONE) and (args.struct_right in NONE)):
 print('Only one structural tag has been set as Non-Existent (NONE). \nCalculation will be stopped. \nPlease set both structural tags as Non-Existent (NONE) if you do not want them in the structure calculation process.')
 sys.exit()

START = time.time()
print('Starting reading input file ...')
if args.variable:
    if args.struct_left in NONE:
        args.struct_left = ''
        LEFT = ''
    else:
        LEFT = args.struct_left
    if args.struct_right in NONE:
        args.struct_right = ''
        RIGHT = ''
    else:
        RIGHT = args.struct_right
    variable_region(filename,args.frequency,args.cutoff)
else:
    LEFT = args.left
    RIGHT = args.right
    if args.major:
        word_major(filename,args.number,args.frequency)
    else:
        word_count(filename,args.number,args.frequency)

CHECK = time.time() - START
print('Finished reading file and fetching aptamers. Processed time: '+str(CHECK)+' sec')
print('Starting secondary structure analysis ...')
os.system('RNAsubopt -e {0} -s -g < subopt_input.txt > subopt_output.txt & wait'.format(args.energy))
CHECK = time.time() - CHECK
print('Finished secondary structure analysis. Processed time: '+str(CHECK)+' sec')
# Create the Dictionary with every secondary structure generated for
# each aptamer. This way, the file subopt_output.txt has no longer to be opened and read again
MFE = 0		# Creates the Minimum Free Energy variable (MFE)
Struct = {}	# Creates the Dictionary for the storage of all the secondary structures
StructScore = {}# Creates the Dictionary for the storage of the scores of each secondary structure
Aptamers = {}	# Creates the Dictionary for the storage of all the Aptamers and their Motifs
HP = {}		# Creates the Dictionary for the Hairpin Scores
BL = {}		# Creates the Dictionary for the Bulge Left Scores
BR = {}		# Creates the Dictionary for the Bulge Right Scores
INT = {}	# Creates the Dictionary for the Intra Strand Scores
QPlex = {}	# Creates the Dictionary for the Quadruplex Scores
Scores = {}	# Creates the Dictionary for the Aptamers Scores
QPL = []        # creates a list to store all the aptamers forming a quadruplex
############### NEW VERSION #############################################
with open('subopt_output.txt','r') as f:
 for line in f:
  if float(line.strip('\n').split(' ')[-1]) < MFE:
   MFE = float(line.strip('\n').split(' ')[-1])
  if line[0].isalpha():
   key = line.strip('\n').split(' ')[0]
   Struct[key] = []
   StructScore[key] = {}
   Aptamers[key] = []
  else:
   dots = str(line.strip('\n').split(' ')[0])
   StructScore[key][dots] = []			# Creates a Dictionary to store motifs for each sequence
   Struct[key].append(line.strip('\n').split(' '))
##########################################################################
# Struct contains now all the structural and energy information of aptamers
# and contains all the aptamers in its keys
# calculate SCORE of ALL Hairpins in the pool
print('Calculating Hairpins motifs scores ...')
hairpin(Struct,MFE,Aptamers,HP,StructScore)
CHECK = time.time() - CHECK
print('Hairpins calculation completed. Processed time: '+str(CHECK)+' sec')
# calculate SCORE of ALL Left Bulges in the pool
print('Calculating Left Bulges motifs scores ...')
left_bulges(Struct,MFE,Aptamers,BL,StructScore)
CHECK = time.time() - CHECK
print('Left Bulges calculation completed. Processed time: '+str(CHECK)+' sec')
# calculate SCORE of ALL Right Bulges in the pool
print('Calculating Right Bulges motifs scores ...')
right_bulges(Struct,MFE,Aptamers,BR,StructScore)
CHECK = time.time() - CHECK
print('Right Bulges calculation completed. Processed time: '+str(CHECK)+' sec')
# calculate SCORE of ALL Intra Strands in the pool
print('Calculating Intra Strands motifs scores ...')
intra(Struct,MFE,Aptamers,INT,StructScore)
CHECK = time.time() - CHECK
print('Intra Strands calculation completed. Processed time: '+str(CHECK)+' sec')
# calculate SCORE of ALL Quadruplex in the pool
print('Calculating Quadruplex motifs scores ...')
quadruplex(Struct,MFE,Aptamers,QPlex,QPL,StructScore)
CHECK = time.time() - CHECK
print('Quadruplex calculation completed. Processed time: '+str(CHECK)+' sec')
# Collates all the single score into a global Dictionary
print('Calculating Aptamers single score ...')
Master_Scores = dict(Counter(HP) + Counter(BL) + Counter(BR) + Counter(INT) + Counter(QPlex))
########## SINGLE APTAMER SCORE CALCULATION #################################
L = list(set(Master_Scores.values()))		#creates a list from all the motif scores
L.sort() 					#and sorts it
Top = percentage(args.upper, len(L))	#calculates the numerosity of the 1% of the motifs list
BEST = []					#sets the list for the Best motifs of the pool
for i in L[-Top:]:				#counts the best 1% motifs
  for key, value in Master_Scores.items():	#then scans of all the motifs
    if value == i:				#and cheks if the motif is in the Best pool
      BEST.append(key)				#and appends to the list

for i in Aptamers.keys():
  if len(Aptamers[i]) == 0:
    Apt_Score = 0
    Scores[i] = [Apt_Score]
  else:
    Match = len(set(Aptamers[i]).intersection(BEST))	#here creates the match list between motifs in the aptamer and those in the best pool
    Apt_Score = (Match/len(Aptamers[i]))*100		#and calculates the score
    Scores[i] = [Apt_Score]				#which is then associated to the aptamer
CHECK = time.time() - CHECK
print('Aptamers Scores calculation completed. Processed time: '+str(CHECK)+' sec')
##################################### PRINTS STRUCTURES SCORES ################
print('Start writing output files ...')
w = csv.writer(open('Structures_Scores.txt', 'w'))
for apt in StructScore.keys():
  for dots in StructScore[apt].keys():
    w.writerow([str(apt)+' '+str(dots)+str(StructScore[apt][dots]).replace("[","").replace("]","").replace("', '"," ").replace("'"," ").replace(",","").replace("  "," ")])
########## PRINTS THE MOTIFS LARGE FILE ########################################
header = 'Aptamer Motif Score'
w = csv.writer(open('Motifs.txt', 'w'))
w.writerow([header])
for key in Aptamers.keys():			
  for value in Aptamers[key]:			# for each motif in each Aptamer
    if value in Master_Scores:			# searches the score of that specific motif
      w.writerow([str(key)+' '+str(value)+' '+str(Master_Scores[value])])
  if key in QPL:
    Scores[key].append('True')
  else:
    Scores[key].append('False')
############ PRINTS APTAMERS SCORE ##############################################
header = 'Aptamer Score Quadruplex'
w = csv.writer(open('Aptamer_Score.txt', 'w'))
w.writerow([header])
for key in Scores:
  w.writerow([str(key)+' '+str(Scores[key][0])+' '+str(Scores[key][1])])
#################################################################################
CHECK = time.time() - START
print('Output file generation completed!')
print('Calculation completed. Thank you for using APTANI². Total Processed time: '+str(CHECK)+' sec')
