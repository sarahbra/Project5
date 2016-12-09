from matplotlib.pyplot import *
C_N = open('crank_nicolson.txt')
T = []
U_CN = np.zeros((21,11))
count = 0

for line in C_N.readlines()[1:-1]:
    temp = line.split(",")
    T.append(float(temp[0]))
    for i in range(len(temp)-1):
        U_CN[count][i] = float(temp[i+1])
    count+=1
C_N.close()

x = np.linspace(0,1,11)
i=0

while i <= count:
    plot(x,U_CN[i,:], label='T= %.4f' %T[i])
    i+=2
legend()
title('Crank-Nicolson diffusion 1dim')
xlabel('x')
ylabel('U(x)')

F_E = open('forward_euler.txt')
U_FE = np.zeros((21,11))
count = 0

for line in F_E.readlines()[1:-1]:
    temp = line.split(",")
    for i in range(len(temp)-1):
        U_FE[count][i] = float(temp[i+1])
    count+=1
F_E.close()

i=0
figure()
while i <= count:
    plot(x,U_FE[i,:], label='T= %.4f' %T[i])
    i+=2
legend()
title('Forward-Euler diffusion 1dim')
xlabel('x')
ylabel('U(x)')

B_E = open('backward_euler.txt')
U_BE = np.zeros((21,11))
count = 0

for line in B_E.readlines()[1:-1]:
    temp = line.split(",")
    for i in range(len(temp)-1):
        U_BE[count][i] = float(temp[i+1])
    count+=1
B_E.close()

i=0
figure()
while i <= count:
    plot(x,U_BE[i,:], label='T= %.4f' %T[i])
    i+=2
legend()
title('Backwards-Euler diffusion 1dim')
xlabel('x')
ylabel('U(x)')

#ana = open('analytic.txt')
#U_ana = np.zeros((21,11))
#count = 0

#for line in ana.readlines()[1:-1]:
#    temp = line.split(",")
#    for i in range(len(temp)-1):
#        U_ana[count][i] = float(temp[i+1])
#    count+=1
#ana.close()

#figure()
#while i <= count:
#    plot(x,U_ana[i,:], label='T= %.2f' %T[i])
#    i+=2

show()
