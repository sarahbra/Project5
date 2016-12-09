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


void printtofile(double t, mat u, int n ){
    ofile << "Time = " << t << endl;
    for (int i=0; i<=n; i++){
        for (int j=0; j<=n; j++) {
            if((j==n)&&(i==n)) {
                ofile << u(i,j);
            } else {
            ofile << u(i,j) << ",";
            }
        }
        ofile << endl;
    }
}

void jakobi_solver (int n, double t_steps, double dx, double dt, mat &A_old) {
    double alpha = dt/(dx*dx);
    mat A_new = zeros<mat>(n+1,n+1);
    for (int t=1;t<=t_steps;t++) {
        for (int i=1;i<n;i++) {
            for (int j=1;j<n;j++) {
                A_new(i,j) = A_old(i,j) + alpha*(A_old(i+1,j)
                    + A_old(i-1,j) - 4*A_old(i,j) + A_old(i,j+1) + A_old(i,j-1));
            }
        }
        A_old = A_new;
        if((t%10==0)||(t==1)) {
            double time = t*dt;
            printtofile(time, A_new,n);
        }
    }
}


double func(double dx, double x_step, double y_step){
    double pi  =3.141592653589793238463;
    return 20*sin(pi*dx*x_step)*sin(pi*dx*y_step);
}

//Found by separation of variables (solving for x while keeping t and y constant and vice versa)
void analytic_Solution (double dx, double dt, double n, double t_steps) {
    mat u = zeros<mat>(t_steps+1,n+1);
    for (int t=1;t<t_steps;t++) {
        for (int i=1;i<n;i++) {
            for (int j=1;j<n;j++) {
                double pi = 3.141592653589793238463;
                u(i,j) = 20*sin(pi*dx*i)*sin(pi*dx*j)*exp(-2*pi*pi*dt*t);
            }
        }
        if((t%10==0)||(t==1)) {
            double time = t*dt;
            printtofile(time, u,n);
        }
    }
}

int main(){
    //Declaring variables
    int n;
    double dx, dt, t_steps;
    char* outfilename;
    outfilename = "test.txt";
    n = 10;
    dx = 0.1;
    dt = dx*dx*0.25;
    t_steps = n*n;
    mat A = zeros<mat>(n+1,n+1);

    for (int i=1;i<n;i++) {
        for (int j=1;j<n;j++) {
            A(i,j) = func(dx,i,j);
        }
    }

    clock_t start, finish;
    start = clock();
    ofile.open(outfilename);
    jakobi_solver(n,t_steps,dx,dt,A);
    ofile.close();
    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;

    outfilename = "analytic.txt";
    ofile.open(outfilename);
    analytic_Solution(dx,dt,n,t_steps);
    return 0;
}
