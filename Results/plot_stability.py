from matplotlib.pyplot import *
C_N1 = open('crank_nicolsondt0.001000.txt')
f = C_N1.readlines()
n = len(f)-2
T1 = []
U_CN1 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T1.append(float(temp[0]))
    if T1[count] > 0.2741 and T1[count] < 0.2759:
        o = count
    if T1[count] > 1.9991 and T1[count] < 2.0009:
        o2 = count
    for i in range(len(temp)-1):
        U_CN1[count][i] = float(temp[i+1])
    count+=1
C_N1.close()

C_N25 = open('crank_nicolsondt0.002500.txt')
f = C_N25.readlines()
n = len(f)-2
T25 = []
U_CN25 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T25.append(float(temp[0]))
    if T25[count] > 0.2741 and T25[count] < 0.2759:
        o_25 = count
    if T25[count] > 1.9991 and T25[count] < 2.0009:
        o2_25 = count
    for i in range(len(temp)-1):
        U_CN25[count][i] = float(temp[i+1])
    count+=1
C_N25.close()

C_N5 = open('crank_nicolsondt0.005000.txt')
f = C_N5.readlines()
n = len(f)-2
T5 = []
U_CN5 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T5.append(float(temp[0]))
    if T5[count] > 0.2741 and T5[count] < 0.2759:
        o_5 = count
    if T5[count] > 1.9991 and T5[count] < 2.0009:
        o2_5 = count
    for i in range(len(temp)-1):
        U_CN5[count][i] = float(temp[i+1])
    count+=1
C_N5.close()

C_N75 = open('crank_nicolsondt0.007500.txt')
f = C_N75.readlines()
n = len(f)-2
T75 = []
U_CN75 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T75.append(float(temp[0]))
    if T75[count] > 0.2741 and T75[count] < 0.2759:
        o_75 = count
    if T75[count] > 1.9991 and T75[count] < 2.0009:
        o2_75 = count
    for i in range(len(temp)-1):
        U_CN75[count][i] = float(temp[i+1])
    count+=1
C_N75.close()

C_N9 = open('crank_nicolsondt0.009000.txt')
f = C_N9.readlines()
n = len(f)-2
T9 = []
U_CN9 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T9.append(float(temp[0]))
    if T9[count] > 0.2741 and T9[count] < 0.2759:
        o_9 = count
    if T9[count] > 1.9991 and T9[count] < 2.0009:
        o2_9 = count
    for i in range(len(temp)-1):
        U_CN9[count][i] = float(temp[i+1])
    count+=1
C_N9.close()

x = np.linspace(0,1,11)
figure()
plot(x,U_CN1[o,:], label='dt=%f' %0.001)
plot(x,U_CN25[o_25,:], label='dt=%f' %0.0025)
plot(x,U_CN5[o_5,:], label='dt=%f' %0.005)
plot(x,U_CN75[o_75,:], label='dt=%f' %0.0075)
plot(x,U_CN9[o_9,:], label='dt=%f' %0.009)
legend()
title('Crank-Nicolson and closed form solution')
xlabel('x')
ylabel('U(x)')
show()
