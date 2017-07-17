#! /usr/bin/python3
import os
with open('suitup_cmd.txt') as f:
	for i in f:
		os.system(i)
