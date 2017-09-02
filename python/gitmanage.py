# -*- coding: utf-8 -*-
# @Author: Helios
# @Date:   2017-08-11 13:22:09
# @Last Modified by:   Helios
# @Last Modified time: 2017-08-14 13:17:44

import re
from subprocess import check_call, check_output

# auto commits to current source without spamming commits
def gitcommit(directory):
	# move to working directory
	check_call(['cd', directory], shell=True)
	# undo last commit if autocommit
	out = check_output(['git', 'log', '-1']).decode('ascii')  #check_output(['git', 'log', '-1'])
	result = re.search(r'(Auto_commit)', out)
	if result is not None:
		# reset latest commit and unstage added files
		check_call(['git', 'reset', 'Head~'])
	# add all unstaged files
	check_call(['git', 'add', '-A'])
	# auto commit
	check_call(['git', 'commit', '-m', 'Auto_commit'])

if __name__ == '__main__':
	import os
	gitcommit(os.getcwd())