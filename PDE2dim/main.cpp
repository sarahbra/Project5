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

double func(double dx, double x_step, double y_step){
    double pi  =3.141592653589793238463;
    return sin(pi*dx*x_step)*sin(pi*dx*y_step);
}

void jakobi_solver (double dx, double dt) {
    double n = 1.0/dx;
    double t_steps = 1000;
    mat A_old = zeros<mat>(n+1,n+1);
    for (int i=1;i<n;i++) {
        for (int j=1;j<n;j++) {
            A_old(i,j) = func(dx,i,j);
        }
    }

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
        if((t%100==0)||(t==1)) {
            double time = t*dt;
            printtofile(time, A_new,n);
        }
    }
}


//Found by separation of variables (solving for x while keeping t and y constant and vice versa)
void analytic_Solution (double dx, double dt) {
    double n = 1.0/dx;
    double t_steps = 1000;
    mat u = zeros<mat>(n+1,n+1);
    for (int t=1;t<t_steps;t++) {
        for (int i=1;i<n;i++) {
            for (int j=1;j<n;j++) {
                double pi = 3.141592653589793238463;
                u(i,j) = 20*sin(pi*dx*i)*sin(pi*dx*j)*exp(-2*pi*pi*dt*t);
            }
        }
        if((t%100==0)||(t==1)) {
            double time = t*dt;
            printtofile(time, u,n);
        }
    }
}

int main(){
    //Declaring variables
    double dx, dt;
    int n, t_steps;
    char* outfilename;
    outfilename = "dx0.1_dt0.6.txt";
    n = 10;
    dx = 0.1;
    dt = dx*dx*0.6;
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
    jakobi_solver(dx,dt);
    ofile.close();
    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;

    outfilename = "ana_dt0.02.txt";
    ofile.open(outfilename);
    analytic_Solution(dx,dt);
    return 0;
}
