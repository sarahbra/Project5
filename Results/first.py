from matplotlib.pyplot import *
C_N1 = open('crank_nicolsondt0.002500.txt')
f = C_N1.readlines()
n = len(f)-2
T1 = []
U_CN1 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T1.append(float(temp[0]))
    if T1[count] > 0.02741 and T1[count] < 0.02759:
        o = count
    if T1[count] > 0.1991 and T1[count] < 0.2009:
        o2 = count
    for i in range(len(temp)-1):
        U_CN1[count][i] = float(temp[i+1])
    count+=1
C_N1.close()
x = np.linspace(0,1,11)



C_N2 = open('backward_eulerdt0.002500.txt')
f = C_N2.readlines()
n = len(f)-2
T2 = []
U_CN2 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T2.append(float(temp[0]))
    if T2[count] > 0.02741 and T2[count] < 0.02759:
        o = count
    if T2[count] > 0.1991 and T2[count] < 0.2009:
        o2 = count
    for i in range(len(temp)-1):
        U_CN2[count][i] = float(temp[i+1])
    count+=1
C_N2.close()

C_N3 = open('forward_eulerdt0.002500.txt')
f = C_N3.readlines()
n = len(f)-2
T3 = []
U_CN3 = np.zeros((n,11))
count = 0

for line in f[1:-1]:
    temp = line.split(",")
    T3.append(float(temp[0]))
    if T3[count] > 0.02741 and T3[count] < 0.02759:
        o = count
    if T3[count] > 0.1991 and T3[count] < 0.2009:
        o2 = count
    for i in range(len(temp)-1):
        U_CN3[count][i] = float(temp[i+1])
    count+=1

ana_x = np.linspace(0,1,100)
temp = 4./np.pi
print T1[0]
t= 0.0025
t2 = 0.65
"""
m_crank_t1 = max(U_CN1[0,:])
m_back_t1 = max(U_CN2[0,:])
m_forw_t1 = max(U_CN3[0,:])
"""
x_ana = np.linspace(0,1,100)
t = 0.
temp = 0
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[0])
"""
m_results_t1 = 1
rel_crank_t1= abs(1-m_crank_t1/m_results_t1)
rel_back_t1 = abs(1-m_back_t1/m_results_t1)
rel_forw_t1 = abs(1-m_forw_t1/m_results_t1)

m_crank_t2 = max(U_CN1[26,:])
m_back_t2 = max(U_CN2[26,:])
m_forw_t2 = max(U_CN3[26,:])

results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[26])
m_results_t2 = max(results)

rel_crank_t2=abs(1-m_crank_t2/m_results_t2)
rel_back_t2 = abs(1-m_back_t2/m_results_t2)
rel_forw_t2 = abs(1-m_forw_t2/m_results_t2)

print"relative error when T = %f :" %T1[0]

print "Forward Euler: %f" %rel_forw_t1
print "Backward Euler: %f" %rel_back_t1
print "Crank-Nicolson: %f" %rel_crank_t1

print"relative error when T = %f :" %T1[26]

print "Forward Euler: %f" %rel_forw_t2
print "Backward Euler: %f" %rel_back_t2
print "Crank-Nicolson: %f" %rel_crank_t2

"""
results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[1])
"""
figure()
plot(x, U_CN1[1,:], label="Crank nicolson")
#plot(ana_x, results, label="Analytical")
title("Cranck-Nicolson and anatlytcal at T=%f" %T1[1])
legend()

results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[26])
figure()
plot(x, U_CN1[26,:], label="Crank nicolson")
#plot(ana_x, results, label="Analytical")
title("Cranck-Nicolson and anatlytcal at T=%f" %T1[26])
legend()

results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[1])

figure()
plot(x, U_CN2[1,:], label="Backward Euler")
#plot(ana_x, results, label="Analytical")
title("Backward Euler and analytical at T=%f" %T1[1])
legend()

results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[26])
figure()
plot(x, U_CN2[26,:], label="Backward Euler")
#plot(ana_x, results, label="Analytical")
title("Backward Euler and analytical at T=%f" %T1[26])
legend()

results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[1])
figure()
plot(x, U_CN3[1,:], label="Forward Euler")
#plot(ana_x, results, label="Analytical")
title("Forward Euler and analytical at T=%f" %T1[1])
legend()

results = (np.sin(np.pi*ana_x))*np.exp(-np.pi**2*T1[26])
figure()
plot(x, U_CN3[26,:], label="Forward Euler")
#plot(ana_x, results, label="Analytical")
title("Forward Euler and analytical at T=%f" %T1[26])
legend()
"""
t = T1[1]
temp = 0
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana , label="T=%f" %t)
figure()
plot(x, U_CN3[1,:], label="Forward Euler")
plot(x_ana, u_ana, label="Analytical")
plot(x, U_CN2[1,:], label="Backward Euler")
plot(x, U_CN1[1,:], label="Crank-Nicolson")
title("Forward Euler, Backward Euler, Crank-Nicolson and analytical at T=%f" %T1[1])
legend()

t = T1[5]
temp = 0
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp

figure()
plot(x, U_CN3[5,:], label="Forward Euler")
plot(x_ana, u_ana, label="Analytical")
plot(x, U_CN2[5,:], label="Backward Euler")
plot(x, U_CN1[5,:], label="Crank-Nicolson")
title("Forward Euler, Backward Euler, Crank-Nicolson and analytical at T=%f" %T1[5])
legend()


figure()
t = 0.
temp = 0
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana , label="T=%f" %t)

temp = 0
t = T1[2]
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana, label="T=%f" %t)

t = T1[5]
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana, label="T=%f" %t)

temp = 0
t = T1[10]
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana, label="T=%f" %t)

temp = 0
t = T1[15]
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana, label="T=%f" %t)

temp = 0
t = T1[18]
for n in range(1,100):
    temp += (2./np.pi)*(-1)**n*(np.sin(n*np.pi*x_ana)/float(n))*np.exp(-(np.pi*n)**2*t)
u_ana = x_ana + temp
plot(x_ana, u_ana, label="T=%f" %t)
title("The time evolution for the analytical solution")
legend()
show()

