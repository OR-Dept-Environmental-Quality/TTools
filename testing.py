from __future__ import division, print_function
from collections import defaultdict

length = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 429.645497802]
azimuths = [x * 360.0 / 8 for x in range(1,8+ 1)]
zone = range(1,int(5))
type = ['LC','ELE']
type2 = type + ['SAMPLE_X','SAMPLE_Y']
type3 = ["SAMPLE_X","SAMPLE_Y"] + type2

#Calculate the unique number of dictionary values
N = len(length) * len(azimuths) * len(zone) * len(type)

# Build the Dictionary
dictkeys = ["LENGTH","SAMPLE_X","SAMPLE_Y","VARIABLE","AZIMUTH","ZONE","VALUE"]
DATA = defaultdict(list)

for d in dictkeys:
	for i in range(0,N):
		DATA[d].append(i)
		
def tree(): return defaultdict(tree)

NODES = tree()

i=0
NODES = tree()
for l in length:
	for a in azimuths:
		for z in zone:
			for t in type2:
				NODES[l][a][z][t] = i
				i = i +1

del(i,x,d,l,t,a,z)

NODE_keys = NODES.keys()
NODE_keys.sort()

####################################################################################################### 
# output shp
# put the dictionary in a list form needed to build the output point feature class
shapekeys = ["LENGTH","SAMPLE_X","SAMPLE_Y","VARIABLE","AZIMUTH","ZONE","VALUE","SAMPLE_X","SAMPLE_Y"]
DATA_shp = [[DATA[k][row] for k in shapekeys] for row in range(0,N)]

AA = [[l,a ,z] for l in length for a in azimuths for z in zone]

NODES_shp = [AA[row] + [NODES[AA[row][0]][AA[row][1]][AA[row][2]][t] for t in type3] for row in range(0,len(AA))]


####################################################################################################### 
# output csv and point files

# Get the stream km dictionary keys and sort them
NODE_keys = NODES.keys()
NODE_keys.sort()

NODES_csv = [[NODES[l][a][z][t] for t in type for a in azimuths for z in zone] for l in length]

# add in the stream km at the beginning of the list
for l in range(0,len(NODE_keys)):
	NODES_csv[l].insert(0,NODE_keys[l])	

# Add the header row
LC_Header = ["Stream_KM"]

for t in range(0,len(type)):
	for a in range(0,len(azimuths)):
		for z in range(0,len(zone)):
			LC_Header.append(type[t]+'_A'+str(a+1)+'_Z'+str(zone[z]))		

NODES_csv.insert(0,LC_Header)
	
i = 0	
