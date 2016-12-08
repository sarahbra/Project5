from matplotlib.pyplot import *
C_N = open('crank_Nicolson.txt')
T = []
U = np.zeros((21,11))
count = 0

for line in C_N.readlines():
    temp = line.split(",")
    T.append(float(temp[0]))
    for i in range(len(temp)-1):
        U[count][i] = float(temp[i+1])
    count+=1
C_N.close()

x = np.linspace(0,1,11)
while i <= count:
    plot(x,U[i,:])
    i+=2
show()
