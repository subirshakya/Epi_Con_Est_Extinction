from __future__ import (division, print_function)
import numpy
from scipy.stats import chisquare as chi

class Node:
	"""This holds the tree data and creates node class"""
	def __init__(self,name="",parent=None,children=None, branchlength = 0):
		self.name = name #Name of node
		self.parent = None #Name of parent (initially set to None, but can be changed)
		if children is None:
			self.children = [] #List to hold children nodes
		else:
			self.children = children #List to append children nodes
		self.brl = branchlength #Branch lengths
		
class Tree:
	"""
	Defines a class of phylogenetic tree, consisting of linked Node objects.
	"""
	
	def __init__(self, data):
		"""
		The constructor really needs to be more flexible, but for now we're 
		going to define the whole tree structure by hand. This just uses
		the same statements we used above. By next Thurs (3/19), see if you can
		write a constructor that takes a parenthetical tree as its argument and 
		builds the corresponding tree in memory. 
		"""

		self.root = Node("root") #Define root
		self.newicksplicer(data, self.root) #Splice newick data
		
	def newicksplicer(self, data, base): 
		""" 
		Splices newick data to create a node based tree. Takes a base argument which is the root node.
		Make sure the file has a single line of newick code. Any miscellaneous characters will throw the program off. Also make sure it does not end in ;"
		Newick should start with ( and end in ). Can modify script to do more.
		Note: The script handles every node as if there is a bifurcation. Any polytomy will not be correctly programmed. 
		"""
		
		data = data.replace(" ", "")[1: len(data)] 	 #Get rid of all spaces and removes first and last parenthesis
		n = 0
		if data.count(",") != 0: #While there is no more comma separated taxa
			for key in range(len(data)): #Find the corresponding comma for a given parenthesis (n will be 0 for the correct comma)
				if data[key] == "(":
					n += 1 #Increase index of n by 1 for 1 step into new node
				elif data[key] == ")":
					n -= 1 #Decrease index of n by 1 for 1 step out node
				elif data[key] == ",": 
					if n == 0: #To check for correct comma
						vals = (data[0:key], data[key+1:len(data)-1]) #Break newick into left and right datasets
						for unit in vals: #For each entry of dataset
							if ":" in unit: #For cases with branch lengths
								data = unit[0:unit.rfind(":")] #get rid of trailing branchlength if provided
								node_creater = Node(data, parent = base) #Create node entry
								node_creater.brl = float(unit[unit.rfind(":")+1:]) #Append branch length of that branch
								base.children.append(node_creater) #Create children
								self.newicksplicer(data, node_creater) #Recursive function
							else: #For case with no branch lengths
								data = unit
								node_creater = Node(data, parent = base)
								base.children.append(node_creater)
								self.newicksplicer(data, node_creater)
						break #Terminate loop, we don't need to look any further
						
	def countsym(self, node):
		"""
		Breaks tree into tow halves at the root node and passes a command to each node.
		"""
		list = [] #List to hold output values
		for child in node.children:
			val = self.countsym2(child)
			nodes = self.nodecount (child)
			list.append((val, nodes))
		return list
		
	def nodecount(self, node):
		"""
		Count number of nodes on tree
		"""
		start = 0
		if node.children == []:
			return 0 #Terminal node returns 0 as there are no nodes aboveit
		else:
			start += 1 #Any non-terminal node will add 1 to the total
			for child in node.children:
				start += self.nodecount(child)
			return start
		
	def length(self, node):
		"""
		Count length to each bifurcation
		"""
		list = []
		if node.children == []: #Terminal branch returns branch length
			return node.brl #Not sure about this part, but without it the numbers are skewed
		else:
			try:
				list.extend(node.brl)
			except:
				list.append(node.brl)
			for child in node.children:
				try:
					list.extend(self.length(child))
				except:
					list.append(self.length(child)) 
			list = filter(lambda x:x != 0, list) #Filter to remove zeros from the list
			return list
		
	def divcount (self):
		"""
		Tallies the lengths to give total diversification rate
		"""
		list = self.length(self.root)
		val = numpy.mean(list)
		return val

	def summarize(self):
		""" 
		Summarize data
		"""
		list = self.countsym(self.root)
		obs = list[0][0]+list[1][0]
		exp = max(list[0][0],list[1][0])*2
		print ("Number of observed species: ", obs)
		print ("Number of extinct species: ", exp-obs)
		print ("% extinct: ", (exp-obs)/(exp)*100)
		print ("chi: ", chi(f_obs = obs, f_exp = exp, ddof = -(obs-1)))

"""		
data = 	#Provide newick here
sim = Tree(data) #Loads tree into program
sim.summarize() #Prints summary stats for the data
"""
