import matplotlib.pyplot as plt
import sys
import pandas
from tqdm import tqdm

def plot(curr_step, DT, filename):
    x = pandas.read_table(filename, header=None, delimiter='\t', usecols=[0])
    y = pandas.read_table(filename, header=None, delimiter='\t', usecols=[1])
    plt.plot(x,y,linewidth=0.75,color='b')
    plt.xlabel('Position', fontsize=12)
    plt.ylabel('Concentration', fontsize=12)
    plt.title("t = " + str(curr_step*DT) + "s")
    plt.ylim(top=1.1, bottom=0)

    filenamePNG = str.replace(filename,'.txt','.png')
    filenamePNG = str.replace(filenamePNG,'data','img')
    plt.savefig(filenamePNG, transparent=False)
    plt.clf()

def getNumber(curr_step):
    if(curr_step < 10):
        return "000" + str(curr_step)
    elif(curr_step < 100):
        return "00" + str(curr_step)
    elif(curr_step < 1000):
        return "0" + str(curr_step)
    else:
        return str(curr_step)

path = ''
DT = 0.0
time_steps = 0

if len(sys.argv) == 4:
    path = sys.argv[1]
    DT = float(sys.argv[2])
    time_steps = int(sys.argv[3])
else:
    print('Error, incorrect call to Python Script')
    exit(1)

filename = 'data/diffusion_1D_'
extension = '.txt'

for i in tqdm(range(time_steps), ncols=100):
    plot(i, DT, filename + getNumber(i) + extension)
