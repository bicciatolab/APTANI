# [Aptamers] deve essere una lista data in input
# [NR] deve essere il numero di strutture da graficare
# [Struct2] è il file da aprire e leggere

import sys,subprocess,os,argparse,csv,re
parser = argparse.ArgumentParser(
    prog='python3.4 grAPhTANI.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=("""

 .----------------. .----------------. .----------------. .----------------. .----------------. .----------------. .----------------. .-----------------..----------------. 
| .--------------. | .--------------. | .--------------. | .--------------. | .--------------. | .--------------. | .--------------. | .--------------. | .--------------. |
| |    ______    | | |  _______     | | |      __      | | |   ______     | | |  ____  ____  | | |  _________   | | |      __      | | | ____  _____  | | |     _____    | |
| |  .' ___  |   | | | |_   __ \    | | |     /  \     | | |  |_   __ \   | | | |_   ||   _| | | | |  _   _  |  | | |     /  \     | | ||_   \|_   _| | | |    |_   _|   | |
| | / .'   \_|   | | |   | |__) |   | | |    / /\ \    | | |    | |__) |  | | |   | |__| |   | | | |_/ | | \_|  | | |    / /\ \    | | |  |   \ | |   | | |      | |     | |
| | | |    ____  | | |   |  __ /    | | |   / ____ \   | | |    |  ___/   | | |   |  __  |   | | |     | |      | | |   / ____ \   | | |  | |\ \| |   | | |      | |     | |
| | \ `.___]  _| | | |  _| |  \ \_  | | | _/ /    \ \_ | | |   _| |_      | | |  _| |  | |_  | | |    _| |_     | | | _/ /    \ \_ | | | _| |_\   |_  | | |     _| |_    | |
| |  `._____.'   | | | |____| |___| | | ||____|  |____|| | |  |_____|     | | | |____||____| | | |   |_____|    | | ||____|  |____|| | ||_____|\____| | | |    |_____|   | |
| |              | | |              | | |              | | |              | | |              | | |              | | |              | | |              | | |              | |
| '--------------' | '--------------' | '--------------' | '--------------' | '--------------' | '--------------' | '--------------' | '--------------' | '--------------' |
 '----------------' '----------------' '----------------' '----------------' '----------------' '----------------' '----------------' '----------------' '----------------' 

An extension tool for the graphic visualization of Aptamer secondary structures generated via APTANI².
grAPhTANI takes as input the 'Structures_Score.csv' file generated via APTANI² in the execution folder. 
Be sure that the file is in the working directory, otherwise grAPhTANI will not work."""))
parser.add_argument("-a","--aptamer", help="Aptamer sequence to be fetched. Can accept more than one",nargs='+')
parser.add_argument("-n","--number",type=int, help="Number of structures to visualize for each aptamer selected. Default is 1", default=1)
parser.add_argument("-t","--algorithm",type=str, help="Visualization algorithm of VARNA. Can choose between radiate, naview, line and circular. Default is naview.", default ="naview")
parser.add_argument("-l","--left",type=str, help="Left tag used during APTANI² calculation. Tags must be submitted in RNA format [AUCG]. Mandatory parameter.", default = "")
parser.add_argument("-r","--right",type=str, help="Left tag used during APTANI² calculation. Tags must be submitted in RNA format [AUCG]. Mandatory parameter.", default = "")
parser.add_argument("-x","--variable",action="store_true",help="Parameter to select if the APTANI² calculation has been performed using the --variable parameter.")
parser.add_argument("-i","--filein",type=str, help="Input file if different from Structure_Scores.csv", default ="")
args = parser.parse_args()
############# VARIABLE DECLARATION #################
yes = ['yes','y', 'ye', '']
no = ['no','n']
Graph = []
Struct2 = []
Aptamers = []
Rankings = {}
img = 1
accepted = ['radiate','line','circular','naview']
############# CHECKS ###############################
if args.algorithm not in accepted:
 print("The chosen algorithm is not supported from VARNA. Please select one algorithm supported between: line, naview, radiate, circular.")
 sys.exit()
if (args.left =="") or (args.right=="") :
 print("No sequences were defined for either left or right tag. Please define your tags according to your experiment.")
 sys.exit()
############# VARIABLE DECLARATION #################
if args.variable:
    LEFT = '^'+args.left
    RIGHT = args.right+'$'
else:
    LEFT = args.left
    RIGHT = args.right
############# APTAMERS FETCHING ####################
for sequence in args.aptamer:
 Aptamers.append(sequence)
############# DATA FETCHING ########################
if args.filein == "":
 with open('Structures_Scores.csv','r') as f:
  for line in f:
   Struct2.append(line.strip('\n').split(' '))
else:
 with open(args.filein,'r') as f:
  for line in f:
   Struct2.append(line.strip('\n').split(' '))
############# DATA ORGANIZATION ####################
for item in Struct2:
 if item[0] in Aptamers:
  if item[0] not in Rankings.keys():
   Rankings[item[0]] = []
#   try:
#    Rankings[item[0]].append(float(item[-1]))
#   except ValueError:
#    continue
#  else:
#   try:
#    if float(item[-1]) not in Rankings[item[0]]:
#     try:
#      Rankings[item[0]].append(float(item[-1]))
#      Rankings[item[0]].sort(reverse=True)
#     except ValueError:
#      continue
#   except ValueError:
#    continue
#  Rankings[item[0]] = Rankings[item[0]][0:args.number]
####################################################
#for item in Struct2:
# try:
#  if item[0] in Rankings.keys():
#   if float(item[-1]) in Rankings[item[0]]:
#    Graph.append(item)
# except ValueError:
#  continue
#for key in Rankings.keys():
# for value in Rankings[key]:
#  for item in Struct2:
#   try:
#    if (item[0] == key) and (float(item[-1]) == value):
#     Graph.append(item)
#     break
#   except ValueError:
#    continue
for key in Rankings.keys():
 for item in Struct2:
  try:
   if item[0] == key:
    Graph.append(item)
  except ValueError:
   continue
###################################################
#choice = input('{0} images will be printed. Continue operation? '.format(len(Graph))).lower()
#if choice in yes:
#   sys.stdout.write('Starting Printing Process')
#elif choice in no:
#   how = input('How many images would like to print? \n')
#   marker = int(how)
####################################################
for test in Graph[:args.number]:
 HP = []
 BG = []
 IS = []
 count = 2
 for i in test[2:-1]:
  count += 1
  if i[0].isalpha() and '|' not in i:				# Hairpin  
   HP.append(test[count])
  elif i[0].isalpha() and (i[-3] == '|') or (i[2] == '|'):	# Bulges
   BG.append(test[count])
   BG.append(test[count+1])
  elif i[0].isalpha() and '|' in i:				# Intras
   IS.append(test[count])
   IS.append(test[count+1])
 if HP:
  HP = '-applyBasesStyle1on "'+(', ').join(HP)+'"'
 else:
  HP = ''
 if BG:
  BG = '-applyBasesStyle2on "'+(', ').join(BG)+'"'
 else:
  BG = ''
 if IS:
  IS = '-applyBasesStyle3on "'+(', ').join(IS)+'"'
 else:
  IS = ''
 title = 'Structure Number: '+str(img)
 if args.variable:
  tags = '-applyBasesStyle4on "'+str(re.search(LEFT,test[0]).start())+'-'+str(re.search(LEFT,test[0]).end())+', '+str(re.search(RIGHT,test[0]).start()+1)+'-'+str(re.search(RIGHT,test[0]).end())+'"'
 else:
  tags = '-applyBasesStyle4on "'+str(re.search(LEFT,test[0]).start()+1)+'-'+str(re.search(LEFT,test[0]).end())+', '+str(re.search(RIGHT,test[0]).start()+1)+'-'+str(re.search(RIGHT,test[0]).end())+'"'
 os.system('java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "{0}" -structureDBN "{1}" -algorithm {2} -autoHelices True -autoInteriorLoops True -backbone "#FF0000" -resolution "3.0" -basesStyle1 "label=#0000e5,fill=#58FAF4,outline=#000000" -basesStyle2 "label=#0000e5,fill=#F7FE2E,outline=#000000" -basesStyle3 "label=#0000e5,fill=#58FA82,outline=#000000" -basesStyle4 "label=#0000e5,fill=#FE2E2E,outline=#000000" {3} {4} {5} {6} -title "{7}" -titleSize 10 -titleColor "#000000" -o {8}'.format(test[0], test[1], args.algorithm, HP, BG, IS, tags, title, 'Image_'+str(img)+':'+str(test[0])+'.png'))
 img+=1
