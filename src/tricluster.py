from statistics import mean
import itertools
class Tricluster:

	def __init__(self, nTimes, nSamples, nPatients):
		self.nTimes = int(nTimes)
		self.nSamples = int(nSamples)
		self.nPatients = int(nPatients)
		self.times = []
		self.samples = []
		self.patients = []

		self.cluster = {}  # format : {(t,s,g): val ...} 

	def addTime(self, time):
		if time not in self.times and len(self.times) < self.nTimes:
			self.times.append(time)
	
	def addSample(self, sample):
		if sample not in self.samples and len(self.samples) < self.nSamples:
			self.samples.append(sample)

	def addPatient(self, patient):
		if len(self.patients) < self.nPatients:
			self.patients.append(patient)

	def addValue(self, time, sample, patient, value):
		if time not in self.times:
			self.times.append(time)
		
		if sample not in self.samples:
			self.samples.append(sample)
		
		if patient not in self.patients:
			self.patients.append(patient)
		
		if (time, sample, patient) not in self.cluster:
			self.cluster[(time, sample, patient)] = value
		
	def hasPatient(self, patient):
		return patient in self.patients 

	def hasTimePoint(self, time):
		return time in self.times
	
	def getTricluster(self):
		return self.cluster

	def getTriclusterCoord(self):
		return self.cluster.keys()
	
	def getPatients(self):
		return self.patients

	def getGPatients(self):
		return list(map(lambda x: int(float(x[2:])), self.patients))

	def getSamples(self):
		return self.samples
	
	def getTimes(self):
		return self.times

	
	def getFeatValues(self, tp, feat):
		vals = filter(lambda x: x[0] == tp and x[1] == feat, self.cluster.keys())

		return [self.cluster[v] for v in vals]

	def getPatientsVals(self, t):
		vals = filter(lambda x: x[0] == t, self.cluster.keys())
		return [self.cluster[v] for v in vals]




	def getSlice(self, g = None, c = None, t = None):
		coords = self.cluster
		if g!= None and c!= None and t!= None:
			return coords[(t,c,g)]
		elif g != None and c != None:
			vals = filter(lambda x: x[2] == g and x[1] == c, coords.keys())
		elif g != None and t != None:
			vals = filter(lambda x: x[2] == g and x[2] == t, coords.keys())
		elif c != None and t != None:
			vals = filter(lambda x: x[1] == c and x[0] == t, coords.keys())
		elif g != None:
			vals = filter(lambda x: x[2] == g, coords.keys())
		elif c != None:
			vals = filter(lambda x: x[1] == c, coords.keys())
		elif t != None:
			vals = filter(lambda x: x[0] == t, coords.keys())
		else:
			return None
		return mygrouper(self.nSamples, [coords[v] for v in vals])

def mygrouper(n, iterable):
	args = [iter(iterable)] * n
	return (list(([e for e in t if e != None] for t in itertools.zip_longest(*args))))


