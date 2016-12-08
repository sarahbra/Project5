#include <ctime>
#include <iostream>
#include <cmath>
#include <fstream>
#include "time.h"
#include<string>
#include<armadillo>
using namespace std;
using namespace arma;

ofstream ofile;

vec GaussElim(double a, vec b, double b_value, double c, int n, vec u, vec v){
    //Forward Substitution
    double m;
    for (int k=2; k<=n; k++) {
        m = a/b(k-1);
        b(k) = b_value - m*c;
        v(k) -= m*v(k-1);
    }

    //Backward Substitution
    u(n)= v(n)/b(n);
    for (int k= n-1; k>0; k--) {
        u(k) = (1.0/b(k))*(v(k) - c*u(k+1));
    }

    u(0) = 0;
    u(n) = 0;
    return u;
}


double func(double dx, double step){
    //double u_i;
    double pi  =3.141592653589793238463;
    return 20*sin(pi*dx*step)*sin(pi*dx*step);

}

void printtofile(vec t, mat u, int n, double t_steps){
    for (int i=0; i<=t_steps; i++){
        ofile << t(i) << ",";
        for (int j=0; j<=n; j++) {
            ofile << u(i,j) << ",";
            if (i==j==n) {
                ofile << u(i,j) << endl;
            }
        }
    }
}

mat backwards_Euler(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    a_value = c_value = -alpha;
    b_value = 1 + 2*alpha;
    mat u_t = zeros<mat>(n+1, n+1);

    vec b(n+1);
    vec u(n+1);
    vec v(n+1);

    //making the vectors
    for (int k=0; k<=n; k++){
        b(k) = b_value;
    }

    for (int k=1; k<n; k++) {
        v(k) = u(k) = func(dx,k);
    }
    //Implementing boundary conditions
    u(n) = v(n) = u(0) = v(0) = 0.0;

    for (int t=1; t<=t_steps; t++) {
        v = GaussElim(a_value, b, b_value, c_value, n, u, v);
        for (int i=1; i<n;i++) {
            u_t(t,i) = v(i);
        }
    }
    return u_t;
}

//Found by separation of variables (solving for x while keeping t and y constant and vice versa)
mat analytic_Solution (double dx, double dt, double n, double t_steps) {
    mat u = zeros<mat>(n+1,n+1);
    for (int t=1;t<=t_steps;t++) {
        for (int i=1;i<n;i++) {
            for (int j=1; j<n;j++) {
                double pi = 3.141592653589793238463;
                u(i,j) = 20*sin(pi*dx*i)*sin(pi*dx*j)*exp(-2*pi*pi*dt*t);
            }
        }
    }
    return u;
}

int main(){
    //Declaring variables
    int n;
    double alpha, dx, dt, t_steps;
    char* outfilename;
    outfilename = "test.txt";
    n = 10;
    dx = 0.1;
    dt = dx*dx*0.25;
    t_steps = n*n;
    vec time(t_steps+1);
    alpha = dt/(dx*dx);

    for (int i=1; i<=t_steps; i++) {
        time(i) = dt*i;
    }

    time(0) = 0.0;

    mat u_xyt = zeros<mat>(n+1, n+1);
    mat u_analytic = zeros<mat>(n+1,n+1);

    ofile.open(outfilename);
    clock_t start, finish;
    start = clock();

    u_xyt = backwards_Euler(n,t_steps,alpha,dx);
    u_xyt += backwards_Euler(n,t_steps,alpha,dx);

    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;

    printtofile(time,u_xyt,n,t_steps);
    ofile.close();

    outfilename = "analytic.txt";
    ofile.open(outfilename);
    u_analytic = analytic_Solution(dx,dt,n,t_steps);
    printtofile(time,u_analytic,n,t_steps);
    ofile.close();
    return 0;
}
