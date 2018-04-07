#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np


s = int(input("Enter the size of matrix\n"))
a_const = 3.0
b_const = 1.0

matr = np.zeros((s+1,s+1))

for i in range(1,s+1):
	matr[i,i] = a_const
for i in range(1,s):
	matr[i,i+1] = (i/s - 0.5) * b_const

for i in range(1,s):
	matr[i+1,i] = matr[i,i+1]

for i in range(1,s+1):
	for k in range(1,s+1):
		if abs(i-k) > 1:
			matr[i,k] = 4/(i+k)**2



matr2 = np.zeros((s,s))
for i in range(s):
	matr2[i,0:s] = matr[i+1,1:s+1] 

#print(matr2)
print()

print(np.linalg.eig(matr2)[0])

