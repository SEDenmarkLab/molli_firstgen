import pandas as pd

def get_xtbout_name(filestr: str):
	"""
	Take xtb output file string and extract name
	"""
	lines = filestr.split('\n')
	name_=''
	for line in lines[60:120]:
		if "coordinate file" in line:
			name_ = line.strip().split()[3].split('.')[0]
			return name_
		else: pass
	raise Exception("Didn't find coordinate file line in output file!")
	
def get_xtbout_sections(filestr):
	"""
	gives tuples with start,stop line indices for sections of output files

	Formatted as dictionary with keys:
	'coeff' for GFN2-xTB coefficients lines
	'wiberg' for wiberg/mayer AO bond orders for each atom
	'fukui' for condensed fukui coef. for each atom
	"""
	lines = filestr.split('\n')
	start_pd = None #Dispersion and polarizability
	end_pd = None
	end_wib = None
	start_wib = None
	fukui_start = None
	fukui_end = None
	count = len(lines) #This is sloppy but will work for now. Fix later.
	for line in lines[::-1]:
		count -= 1 #Corrects for zero indexing, and only has to be performed at the top of the for loop, not after each condition
		idx = count
		if 'Mol. C6AA' in line:
			end_pd = idx-2
		elif 'covCN' in line:
			start_pd = idx+1
		elif 'Topologies differ' in line:
			end_wib = idx-3
		elif 'Z sym  total' in line:
			start_wib = idx+2
		elif 'f(+)' in line:
			fukui_start = idx+1
			break
		elif 'Property Printout' in line:
			fukui_end = idx-2
		elif '(HOMO)' in line:
			homoline = line
		elif '(LUMO)' in line:
			lumoline = line
		else: pass
	outdict = {
		'coeff':(start_pd,end_pd),
		'wiberg':(start_wib,end_wib),
		'fukui':(fukui_start,fukui_end),
	}
	checklist = [start_pd,start_wib,fukui_start,end_pd,end_wib,fukui_end]
	if None in checklist:
		print(checklist)
		raise Exception('Did not find one of the bounds in file!')
	return outdict

def get_xtb_coef(filestr: str,outdict: dict):
	"""
	This is meant to operate with dictionary of file sections from 
	get_xtbout_sections()
	"""
	lines = filestr.split('\n')
	st,end = outdict['coeff']
	end+=1
	# print(outdict['coeff'])
	# print(st,end,type(st),type(end))
	# print(lines)
	section = lines[st:end]
	outdict = {}
	for line in section:
		split = line.strip().split()
		outdict[int(split[0])] = {'symbol':split[2],'disp':split[5],'pol':split[6],'charge':split[4],'covCN':split[3]}
	df = pd.DataFrame.from_dict(outdict,orient='index')
	return df

def get_xtb_fukui(filestr: str,outdict: dict):
	"""
	This pulls out fukui incices and writes them to a dataframe
	"""
	lines = filestr.split('\n')
	st,end = outdict['fukui']
	end+=1
	# print(outdict['coeff'])
	# print(st,end,type(st),type(end))
	# print(lines)
	section = lines[st:end]
	outdict = {}
	for line in section:
		split = line.strip().split()
		atomid = int(''.join([f for f in split[0] if f.isdigit()]))
		outdict[atomid] = {'f+':split[1],'f-':split[2],'f0':split[3]}
	df = pd.DataFrame.from_dict(outdict,orient='index')
	return df

def get_xtb_wiberg(filestr: str,outdict: dict):
	"""
	Get MAX wiberg bond order for each atom
	"""
	lines = filestr.split('\n')
	st,end = outdict['wiberg']
	end+=1
	# print(outdict['coeff'])
	# print(st,end,type(st),type(end))
	# print(lines)
	section = lines[st:end]
	outdict = {}
	for line in section:
		split = line.strip().split()
		if len(split) == 3 or len(split) == 6: continue
		else:
			outdict[int(split[0])] = {'max_bond_order':split[7]}
	df = pd.DataFrame.from_dict(outdict,orient='index')
	return df

def extract_xtb_atomic_properties(xtbout: str,xtb_coef = True,fukui = True,wiberg=True):
	"""
	Pass in output files from xtb optimizations and parse:

	charges
	largest wiberg bond orders for each atom
	polarizabilitis 
	fukui indices
	"""
	# with open(xtbout,'rb') as g:
	# 	filestr = g.read().decode('utf-8')
	filestr = xtbout
	sections = get_xtbout_sections(filestr)
	out = []
	if xtb_coef == True: 
		df1 = get_xtb_coef(filestr,sections)
		out.append(df1)
	if fukui == True: 
		df2 = get_xtb_fukui(filestr,sections)
		out.append(df2)
	if wiberg == True: 
		df3 = get_xtb_wiberg(filestr,sections)
		out.append(df3)
	outdf = pd.concat(out,axis=1)
	name = get_xtbout_name(filestr)
	outdf.name = name
	return outdf
