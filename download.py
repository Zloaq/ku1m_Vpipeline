#!/Users/motomo/opt/anaconda3/bin/python3

import pexpect
import os
import sys


def download(command, password):

	child = pexpect.spawn(command)

	while True:
		i = child.expect(['password:', pexpect.EOF, pexpect.TIMEOUT])
		if i == 0:
			child.sendline(password)
		elif i == 1:
			break
		elif i == 2:
			print("SCP command timed out.")
			break

	print(child.before.decode())



def sort(obnamelist, obname):
	
	script_dir = os.path.dirname(os.path.abspath(__file__))
	file_path = os.path.join(script_dir, 'lib2', obname)
	if not os.access(file_path, os.R_OK):
		print(file_path, ' is not found')
		sys.exit()
	with open(file_path, 'r') as f1:
		for line in f1:
			if line.startswith('search_names'):
				varr = line.split()
				break
	if len(fitslist) != len(obnamelist):
		print('list is irregal')
		sys.exit()
	
	for i1 in varr[1:]:
		index1 = []
		index2 = []
		for i2 in range(len(obnamelist)):
			if i1 in obnamelist[i2]:
				index1.append(i2)
				print(obnamelist[i2], ' key---', i1)
			else:
				index2.append(i2)
		if len(index1) != 0:
			break
	
	return index1, index2
