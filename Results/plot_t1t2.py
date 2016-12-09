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

ana = open('analytic.txt')
U_ana = np.zeros((21,11))
count = 0

for line in ana.readlines():
    temp = line.split(",")
    for i in range(len(temp)-1):
        U_ana[count][i] = 4.0/np.pi*float(temp[i+1])
    count+=1
ana.close()

figure()
plot(x,U_CN[8,:], label='CN T=%.4f' %T[2])
plot(x,U_ana[8,:], label='Analytic T=%.4f' %T[2])
legend()
title('Crank-Nicolson and closed form solution')
xlabel('x')
ylabel('U(x)')

figure()
plot(x,U_CN[18,:], label='CN T=%.4f' %T[18])
plot(x,U_ana[18,:], label='Analytic T=%.4f' %T[18])
legend()
title('Crank-Nicolson and closed form solution')
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

figure()
plot(x,U_FE[8,:], label='FE T=%.4f' %T[4])
plot(x,U_ana[8,:], label='Analytic T=%.4f' %T[4])
legend()
title('Forward-Euler and closed form solution')
xlabel('x')
ylabel('U(x)')

figure()
plot(x,U_FE[18,:], label='FE T=%.4f' %T[18])
plot(x,U_ana[18,:], label='Analytic T=%.4f' %T[18])
legend()
title('Forward-Euler and closed form solution')
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

figure()
plot(x,U_BE[8,:], label='BE T=%.4f' %T[4])
plot(x,U_ana[8,:], label='Analytic T=%.4f' %T[4])
legend()
title('Backwards-Euler and closed form solution')
xlabel('x')
ylabel('U(x)')

figure()
plot(x,U_BE[18,:], label='BE T=%.4f' %T[18])
plot(x,U_ana[18,:], label='Analytic T=%.4f' %T[18])
legend()
title('Backwards-Euler and closed form solution')
xlabel('x')
ylabel('U(x)')

show()
