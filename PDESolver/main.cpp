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

void initializePrint(char* method, int n,double dx,double dt, int time_steps) {
    char outfilename[60];
    sprintf(outfilename, "%sdt%f.txt", method, dt);


    //outfilename = method;
    ofile.open(outfilename);
    ofile << method << "   " << "n= " << n << " dx= " << dx << " dt= " << dt << "t_steps " << time_steps << endl;
}

void finalizePrint(double seconds){
    ofile << "Computation took s= " << seconds << " seconds" << endl;
    ofile.close();
}

void printtofile(double t, vec u, int n ){
    ofile << t << ",";
    for (int i=0; i<n; i++){
        ofile << u(i) << ",";
    }
    ofile << u(n) << endl;
}


vec GaussElim(double a, vec b, double b_value, double c, int n, vec u, vec v){
    //Forward Substitution
    double m;
    cout << "n: " << n << endl;
    for (int k=2; k<=n; k++) {
        cout << "k " << k << endl;
        m = a/b(k-1);
        b(k) = b(k) - m*c;
        v(k) = v(k) - m*v(k-1);
    }

    //Backward Substitution
    u(n)= v(n)/b(n);
    cout << "u(n)" << u(n) << endl;
    for (int k= n-1; k>0; k--) {
        u(k) = (1.0/b(k))*(v(k) - c*u(k+1));
        cout << "u(k)" << u(k) <<"   k " << k << endl;

    }

    u(0) = 0;
    u(n) = 1;
    return u;
}


double func(double dx, double step){
    //double u_i;
    double pi  =3.141592653589793238463;
    //return sin(pi*dx*step);
    return 0;

}

vec forward_step(double n, double alpha, vec u, vec unew) {
    for (int i=1; i<n; i++) {
        unew(i) = alpha*u(i-1) + (1-2*alpha) * u(i) + alpha*u(i+1);
    }
    return unew;
}

void forward_Euler(int n, int t_steps, double alpha, double dx, double dt) {
    char* method = "forward_euler";

    initializePrint(method, n, dx, dt, t_steps);
    vec u(n+1);
    vec unew(n+1);

    u(0) = unew(0) = 0.0;//
    u(n) = unew(n) = 1.0;

    //making the vectors
    for (int k=1; k<n; k++){
        u(k) = func(dx, k);
        unew(k) = 0;
    }

    for (int t=1;t<=t_steps;t++) {
        u = forward_step(n, alpha, u, unew);
        if(t%10==0 || t==1) {
            double time = t*dt;
            printtofile(time, u, n);
        }
    }
}

void backwards_Euler(int n, int t_steps, double alpha, double dx, double dt) {
    char* method = "backward_euler";
    initializePrint(method, n, dx, dt, t_steps);
    double a_value, c_value, b_value;
    a_value = c_value = -alpha;
    b_value = 1 + 2*alpha;

    vec b(n+1);
    vec u(n+1);
    vec v(n+1);

    //making the vectors
    for (int k=0; k<=n; k++){
        b(k) = b_value;
    }

    for (int k=1; k<n; k++) {
        u(k) = v(k) =func(dx,k);
    }
    //Implementing boundary conditions
    u(n) = v(n) = 1.0;
    u(0) = v(0) = 0.0;

    for (int t=1; t<=t_steps; t++) {
        v = GaussElim(a_value, b, b_value, c_value, n, u, v);
        if (t%10==0 || t==1) {
            double time = t*dt;
            printtofile (time,v,n);
        }

    }
}

void crank_Nicolson(int n, int t_steps, double alpha, double dx, double dt) {
    char* method = "crank_nicolson";
    initializePrint(method, n, dx, dt, t_steps);
    double a_value, c_value, b_value;
    vec b(n+1);
    a_value = c_value = -alpha;
    b_value = 1 + 2*alpha;

    vec u(n+1);
    vec v(n+1);

    //making the vectors
    for (int k=1; k<n; k++){
        b(k) = b_value;
        u(k) = func(dx,k);
    }

    b(0) = b(n) = b_value;
    u(0) = 0.0;
    u(n) = 1.0;

    //GaussElim(a_value,b, b_value,c_value,u,n,v);
    for (int t=1;t<=t_steps;t++) {
        v = forward_step(n,alpha,u,v);
        v(0) = 0.0;
        v(n) = 1.0;
        u = GaussElim(a_value,b, b_value,c_value,n,u,v);
        if (t%10==0 || t==1) {
            double time = t*dt;
            printtofile(time,v,n);
        }
    }
}

void analytic_Solution (double dx, double dt, double n, double t_step) {
    vec u(n+1);
    u(0) = u(n) = 0.0;
    for (int t=1;t<=t_step;t++) {
        for (int i=1;i<n;i++) {
            double pi = 3.141592653589793238463;
            u(i) = sin(pi*dx*i)*exp(-pi*pi*dt*t);
        }
        if (t%10==0 || t==1) {
            double time = t*dt;
            printtofile(time,u,n);
        }
    }
}

int main(){
    //Declaring variables
    int n;
    double alpha, dx, t_steps, dt, final_t, t, seconds;
    vec dt_list(5);
    n = 10;
    dx = 0.1;
    //dt = dx*dx*0.25;
    //t_steps = 2*n*n;
    dt_list(0) = dx*dx*0.01;
    dt_list(1) = dx*dx*0.1;
    dt_list(2) = dx*dx*0.25;
    dt_list(3) = dx*dx*0.35;
    dt_list(4) = dx*dx*0.5;
    final_t = 1;

    clock_t start, finish;

    for (int i = 0; i<5 ; i++){
        t_steps = final_t/dt_list(i);

        start = clock();
        dt = dt_list(i);
        alpha = dt/(dx*dx);
        forward_Euler(n,t_steps,alpha,dx ,dt);

        finish =clock();
        double t = ((finish-start));
        double seconds = t/CLOCKS_PER_SEC;
        finalizePrint(seconds);
    }

    for (int i = 0; i<5; i++){
        t_steps = final_t/dt_list(i);
        dt = dt_list(i);
        alpha = dt/(dx*dx);
        start = clock();
        backwards_Euler(n,t_steps,alpha,dx, dt);
        finish = clock();
        t = ((finish-start));
        seconds = t/CLOCKS_PER_SEC;
        finalizePrint(seconds);
    }

    for(int i = 0; i<5; i++){
        t_steps = final_t/dt_list(i);
        dt = dt_list(i);
        alpha = dt/(dx*dx);
        start = clock();
        crank_Nicolson(n,t_steps,alpha/2,dx, dt);

        finish =clock();
        t = ((finish-start));
        seconds = t/CLOCKS_PER_SEC;

        finalizePrint(seconds);
    }
    dt = 0.1*0.1*0.25;
    char* outfilename_ana = "analytic.txt";
    ofile.open(outfilename_ana);
    analytic_Solution(dx,dt,n,t_steps);
    ofile.close();
    return 0;
}
