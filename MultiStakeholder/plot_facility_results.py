# Facility Location Case Study
# Created by Alex Dowling (adowling2@wisc.edu)
# University of Wisconsin-Madison
# Department of Chemical and Biological Engineering
# Last Modified Feb. 17th, 2016

# This Python script important the results from text files and plots the solution

### Import

import re

import matplotlib.pyplot as plt
import matplotlib.text as mpl_text
import matplotlib.font_manager as fm

### Legend Definitions

class AnyObject(object):
    def __init__(self, text, color, size):
        self.my_text = text
        self.my_color = color
        self.my_size = size

class AnyObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#        print orig_handle
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        myfont = fm.FontProperties(fname="/usr/share/fonts/truetype/ttf-ancient-scripts/Symbola605.ttf")
        patch = mpl_text.Text(x=0, y=0, text=orig_handle.my_text, color=orig_handle.my_color, verticalalignment=u'baseline', 
                                horizontalalignment=u'left', multialignment=None, 
                                rotation=0, linespacing=None, 
                                rotation_mode=None,fontproperties=myfont,
                                size = orig_handle.my_size)
        handlebox.add_artist(patch)
        return patch

### Data type definitions

class ProbData(object):
	def __init__(self, farmLoc, cityLoc, waterLoc, candLoc, xCand, yCand):
		self.farmLoc = farmLoc
		self.cityLoc = cityLoc
		self.waterLoc = waterLoc
		self.candLoc = candLoc
		self.xCand = xCand
		self.yCand = yCand

class ResultsData(object):
	def __init__(self, smFac, lrgFac, ffm, basename):
		self.smFac = smFac
		self.lrgFac = lrgFac
		self.ffm = ffm
		self.basename = basename

### Function definitions

def readDataFile(fileName):

	fileD = open(fileName,'r')

	# Parse farms first
	ln = fileD.readline()
	r = re.match(r'Farms: ([0-9]+)', ln)
	nFarms = int(r.group(1))

	farmLoc = []

	for i in range(nFarms):
		ln = fileD.readline()
		r = re.match(r'([0-9\.]+),([0-9\.]+)', ln)
		farmLoc.append([float(r.group(1)), float(r.group(2))])

	# Skip empty line
	ln = fileD.readline()
	
	# Parse cities next
	ln = fileD.readline()
	r = re.match(r'Cities: ([0-9]+)', ln)
	nUrban = int(r.group(1))

	cityLoc = []

	for i in range(nUrban):
		ln = fileD.readline()
		r = re.match(r'([0-9\.]+),([0-9\.]+)', ln)
		cityLoc.append([float(r.group(1)), float(r.group(2))])

	# Skip empty line
	ln = fileD.readline()
	
	# Parse water next
	ln = fileD.readline()
	r = re.match(r'Water: ([0-9]+)', ln)
	nWater = int(r.group(1))

	waterLoc = []

	for i in range(nWater):
		ln = fileD.readline()
		r = re.match(r'([0-9\.]+),([0-9\.]+)', ln)
		waterLoc.append([float(r.group(1)), float(r.group(2))])

	# Skip empty line
	ln = fileD.readline()

	# Parse facilities next
	ln = fileD.readline()
	r = re.match(r'Facilities: ([0-9]+)', ln)
	nCand= int(r.group(1))

	candLoc = []

	for i in range(nCand):
		ln = fileD.readline()
		r = re.match(r'([0-9\.]+),([0-9\.]+)', ln)
		candLoc.append([float(r.group(1)), float(r.group(2))])

	# Skip empty line
	ln = fileD.readline()

	fileD.close()

	# Copy candidate locations into two lists suitable for scatter plot construction
	xCand = []
	yCand = []
	for i in range(nCand):
		xCand.append(candLoc[i][0])
		yCand.append(candLoc[i][1])

	return ProbData(farmLoc,cityLoc,waterLoc,candLoc, xCand, yCand)


def readResultsFile(fileName):

	fileR =  open(fileName,'r')
	ln = fileR.readline()

	smFac = []
	lrgFac = []
	ffm = []

	mode = 0

	while ln:
		if ln == 'Small Facilities:\n':
			mode = 1
		elif ln == 'Large Facilities:\n':
			mode = 2
		elif ln == 'Farm-Facility Matrix (Sparse Form):\n':
			mode = 3 
		elif mode == 1:
			smFac.append(int(ln))
		elif mode == 2:
			lrgFac.append(int(ln))
		elif mode == 3:
			if ln[0:3] == 'z =':
				print ln
			elif ln[0:3] == 'f =':
				print ln
			elif ln[0:7] == 'alpha =':
				print ln
			elif ln[0:6] == 'CVaR =':
				print ln
			else:
				r = re.match(r'\(([0-9]+),([0-9]+),([0-9\.]+)\)', ln)
				ffm.append([int(r.group(1)),int(r.group(2)),float(r.group(3))])
			
		ln = fileR.readline()
	
	fileR.close()
	
	# Parse base file name - requires no "." in file name except before extension
	r = re.match(r'([^\.]+)\.txt',fileName)
	
	return ResultsData(smFac, lrgFac, ffm,r.group(1))


def plotResults(pd, rd, markCandLoc):

	### Setup
	
	# Specify unicode characters for plot symbols
	cow = u'\U0001F404'
	#city = u'\U0001F3E2'
	#city = u'\U0001F3D8'
	city = u'\U0001F3E0'
	water = u'\U0001F4A7'
	sFac = u'\u2699'
	lFac = u'\U0001F3ED'

	# Specify standard, small and large symbol sizes
	stnd = 18
	sm = 14
	lrg = 22

	### Initialization

	# Create figure
	fig = plt.figure()
	ax = fig.add_subplot(111)

	### Plot Problem Data

	# Plot facility candidate locations
	cl = ax.scatter(x=pd.xCand,y=pd.yCand,marker=".",color="grey",label="Candidate Locations")

	# Draw unit box
	ax.plot([0,0],[0,1],color="black")
	ax.plot([0,1],[1,1],color="black")
	ax.plot([1,1],[1,0],color="black")
	ax.plot([1,0],[0,0],color="black")
	
	# Adjust axes
	plt.xlim(-0.05,1.05)
	plt.ylim(-0.05,1.05)

	# Specify font for unicode characters
	myfont = fm.FontProperties(fname="/usr/share/fonts/truetype/ttf-ancient-scripts/Symbola605.ttf")

	# Plot farms
	for i in range(len(pd.farmLoc)):
		ax.text(pd.farmLoc[i][0],pd.farmLoc[i][1],cow,color="red",va="center",ha="center",fontproperties=myfont,fontsize=stnd)

	# Plot urban centers
	for i in range(len(pd.cityLoc)):
		ax.text(pd.cityLoc[i][0],pd.cityLoc[i][1],city,color="green",va="center",ha="center",fontproperties=myfont,fontsize=stnd)

	# Plot water sheds
	for i in range(len(pd.waterLoc)):
		ax.text(pd.waterLoc[i][0],pd.waterLoc[i][1],water,color="blue",va="center",ha="center",fontproperties=myfont,fontsize=stnd)

	### Plot results

	# Draw small facilities
	for i in range(len(rd.smFac)):
		ax.text(pd.candLoc[rd.smFac[i]-1][0],pd.candLoc[rd.smFac[i]-1][1],sFac,color="black",va="center",ha="center",fontproperties=myfont,fontsize=sm)

	# Draw large facilities
	for i in range(len(rd.lrgFac)):
		ax.text(pd.candLoc[rd.lrgFac[i]-1][0],pd.candLoc[rd.lrgFac[i]-1][1],lFac,color="black",va="center",ha="center",fontproperties=myfont,fontsize=lrg)

	# Draw non-zero elements of farm-to-facility matrix
	for i in range(len(rd.ffm)):
		ax.plot([pd.farmLoc[rd.ffm[i][0] - 1][0], pd.candLoc[rd.ffm[i][1] - 1][0]], [pd.farmLoc[rd.ffm[i][0] - 1][1], pd.candLoc[rd.ffm[i][1] - 1][1]], color = "grey", linestyle='--')

	# Mark certain candidate locations
	for i in markCandLoc:
		circ = plt.Circle((pd.candLoc[i - 1][0],pd.candLoc[i - 1][1]),.05,color='m',fill=False)
		fig.gca().add_artist(circ)

	### Create legend

	obj_0 = AnyObject(cow, "red",stnd)
	obj_1 = AnyObject(city, "green",stnd)
	obj_2 = AnyObject(water, "blue",stnd)
	obj_3 = AnyObject(sFac, "black",sm)
	obj_4 = AnyObject(lFac, "black",lrg)


	plt.legend([obj_0, obj_2, obj_1, cl, obj_3, obj_4], ['Farms', 'Water Sheds', 'Urban Centers', 'Candidate Locations', 'Small Facilities','Large Facilities'],
           handler_map={obj_0:AnyObjectHandler(), obj_1:AnyObjectHandler(), obj_2:AnyObjectHandler(), obj_3:AnyObjectHandler(), obj_4:AnyObjectHandler()},
           bbox_to_anchor=(0.5, 1.125), loc="upper center", ncol=3, fancybox="True", shadow="True")

	### Save plots
	
	fig.savefig(rd.basename + '.eps')
	fig.savefig(rd.basename + '.png')

#########################
### Script part

pd = readDataFile('test_data.txt')

resultsFileNames = []


iSaveAlt = [1, 51, 73, 88, 96, 99]
iSaveTrad = [1, 61, 65, 76, 91, 92, 94]

iStakeholders = [41, 25, 35, 45]

# Tip. In order to calculate alpha from filename...
# alpha = (iSaveTrad or iSaveAlt - 1)/100 

for i in range(4):
	for j in range(2):
		if j == 0:
			t = 'alt'
		else:
			t = 'trad'
		resultsFileNames.append('Solution_Obj_' + str(i+1) + '_only_' + t + '.txt')
		
for j in range(len(iSaveTrad)):
	resultsFileNames.append('Solution_CVaR_' + str(iSaveTrad[j]) + '_Trad_Nadir.txt')
	
for j in range(len(iSaveAlt)):
	resultsFileNames.append('Solution_CVaR_' + str(iSaveAlt[j]) + '_Alt_Nadir.txt')

for j in range(len(iStakeholders)):
	resultsFileNames.append('Solution_Stakeholder_' + str(iStakeholders[j]) + '_Alt_Nadir.txt')

results = []
	
for i in range(len(resultsFileNames)):
	rd = readResultsFile(resultsFileNames[i])
	results.append(rd)

lrgFacLoc_trad = []
lrgFacLoc_alt = []

i1 = 8
i2 = i1 + len(iSaveTrad)
i3 = i2 + len(iSaveAlt)
i4 = i3 + len(iStakeholders)


# Determine unique large facility locations with traditional nadir point
for i in range(i1,i2):
	for j in results[i].lrgFac:
		if j not in lrgFacLoc_trad:
			lrgFacLoc_trad.append(j)

# Determine unique large facility locations with alternate nadir point				
for i in range(i2,i3):
	for j in results[i].lrgFac:
		if j not in lrgFacLoc_alt:
			lrgFacLoc_alt.append(j)
				
for i in range(i1):
	plotResults(pd,results[i],[])

for i in range(i1,i2):
	plotResults(pd,results[i],lrgFacLoc_trad)
	
for i in range(i2,i3):
	plotResults(pd,results[i],lrgFacLoc_alt)

for i in range(i3,i4):
	plotResults(pd,results[i],lrgFacLoc_alt)

plt.show()


