from collections import defaultdict

length = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 429.645497802]
azimuths = [x * 360.0 / 8 for x in range(1,8+ 1)]
zone = range(0,int(5))
type = ['LC','ELE']

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

type2 = type + ['SAMPLE_X','SAMPLE_Y']

i=0
NODES = tree()
for l in length:
	for a in azimuths:
		for z in zone:
			for t in type2:
				NODES[l][a][z][t] = i
				i = i +1

del(i,x,d,l,t,a,z)

from operator import itemgetter
NODE_keys = NODES.keys()
NODE_keys.sort()


# put the dictionary in a list form needed to build the output point feature class
shapekeys = ["LENGTH","SAMPLE_X","SAMPLE_Y","VARIABLE","AZIMUTH","ZONE","VALUE","SAMPLE_X","SAMPLE_Y"]
DATA_shp = [[DATA[k][row] for k in shapekeys] for row in range(0,N)]
#NODES_shp = [[NODES[l][t][a][z]] for l in length for t in type for a in azimuths for z in zone]

#NODES_shp = []
#for i in range(0,N):
	#for l in length:
		#for t in type:
			#for a in azimuths:
				#for z in zone:
					#NODES_shp[i].append(NODES[l][t][a][z])


# put into format for csv
NODES_csv = [[NODES[l][a][z][t] for t in type for a in azimuths for z in zone] for l in length]

del(l,t,a,z)
for l in range(0,len(NODE_keys)):
	NODES_csv[l].insert(0,NODE_keys[l])

lcDataColHeaders = ["Stream_KM"]

# Add the header row
for t in range(0,len(type)):
	for a in range(0,len(azimuths)):
		for z in range(0,len(zone)):
			lcDataColHeaders.append(type[t]+'_A'+str(a+1)+'_Z'+str(zone[z]))		

NODES_csv.insert(0,lcDataColHeaders)
	
i = o	