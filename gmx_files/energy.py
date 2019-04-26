#! /home/aware/anaconda3/bin/ipython

# coding: utf-8
# by Artemii Yanushevskyi

import os, sys
import pandas as pd
import matplotlib.pyplot as plt



try:
    case = sys.argv[1]
    f = open(os.path.join(case, 'energy.xvg'))
except FileNotFoundError as ex:
    print("File 'energy.xvg' doesn't exist. We need to convert the binary energy file (ener.edr) to energy text file (energy.xvg) using command `gmx energy` and typing '3 4 5 6'.")
    raise ex
except IndexError as ex:
    print("Write the version of the simulation.")
    raise ex

text = f.read()

X = []

# parsing
for line in text.split('\n'):
    if line != '':
        if line[0] in '#@' or line == '':
            if line[0:3] == '@ s':
                print(line[4:])
            continue

        values = map(float, line.split())
        print(line)
        X.append(list(values))


df = pd.DataFrame(X, columns=['index', 'Potential', 'Kinetic En.', 'Total Energy', 'Temperature']).drop('index', axis=1)
df.plot(ls='-', marker='.')

location = os.path.join(case, 'energies.png')
plt.savefig(location)
os.system("xdg-open " + location)
