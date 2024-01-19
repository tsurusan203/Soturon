import matplotlib.pyplot as plt

x_list=[]

f=open('origin.dat', 'rt')
for line in f:
    data=line[:-1].split('')
    x_list.append(float(data[0]))

plt.plot(x_list, 'o',linestyle='None')

plt.xlabel('x')

plt.grid(True)

plt.show()