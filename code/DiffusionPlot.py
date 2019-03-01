import matplotlib.pyplot as plt
import sys
import pandas

filename = ''
time = 0

if len(sys.argv) == 3:
    filename = sys.argv[1]
    time = sys.argv[2]
else:
    print('Error, incorrect call to Python Script')
    exit(1)

x = pandas.read_table(filename, header=None, delimiter='\t', usecols=[0])
y = pandas.read_table(filename, header=None, delimiter='\t', usecols=[1])
plt.plot(x,y,linewidth=0.75,color='b')
plt.xlabel('Position', fontsize=12)
plt.ylabel('Concentration', fontsize=12)
plt.title("t = " + str(time))
plt.ylim(top=1.1, bottom=0)

filenamePNG = str.replace(filename,'.txt','.png')
filenamePNG = str.replace(filenamePNG,'data','img')
plt.savefig(filenamePNG, transparent=False)