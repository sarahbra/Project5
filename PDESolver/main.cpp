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

void initializePrint(string method, int n,double dx,double dt) {
    string outfilename;
    outfilename = method + ".txt";
    ofile.open(outfilename);
    ofile << method << "   " << "n= " << n << " dx= " << dx << " dt= " << dt << endl;
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
    return 20*sin(pi*dx*step);

}




vec forward_step(double n, double alpha, vec u, vec unew) {
    for (int i=1; i<n; i++) {
        unew(i) = alpha*u(i-1) + (1-2*alpha) * u(i) + alpha*u(i+1);
    }
    return unew;
}

void forward_Euler(int n, int t_steps, double alpha, double dx, double dt) {
    string method = "forward_euler";
    initializePrint(method, n, dx, dt);
    vec u(n+1);
    vec unew(n+1);

    u(0) = unew(0) = u(n) = unew(n) = 0.0;

    //making the vectors
    for (int k=1; k<n; k++){
        u(k) = func(dx, k);
        unew(k) = 0;
    }

    for (int t=1;t<=t_steps;t++) {
        u = forward_step(n, alpha, u, unew);
        if(t%10==0 || t==1) {
            double time = t*dx*dx*0.25;
            printtofile(time, u, n);
        }
    }
}

void backwards_Euler(int n, int t_steps, double alpha, double dx, double dt) {
    string method = "backward_euler";
    initializePrint(method, n, dx, dt);
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
    u(n) = v(n) = u(0) = v(0) = 0.0;


    for (int t=1; t<=t_steps; t++) {
        v = GaussElim(a_value, b, b_value, c_value, n, u, v);
        if (t%10==0 || t==1) {
            double time = t*dx*dx*0.25;
            printtofile (time,v,n);
        }

    }
}

void crank_Nicolson(int n, int t_steps, double alpha, double dx, double dt) {
    string method = "crank_nicolson";
    initializePrint(method, n, dx, dt);
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
    u(0) = u(n) = 0.0;

    //GaussElim(a_value,b, b_value,c_value,u,n,v);
    for (int t=1;t<=t_steps;t++) {
        v = forward_step(n,alpha,u,v);
        v(0) = v(n) = 0;
        u = GaussElim(a_value,b, b_value,c_value,n,u,v);
        if (t%10==0 || t==1) {
            double time = t*dx*dx*0.25;
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
            u(i) = 20*sin(pi*dx*i)*exp(-pi*pi*dt*t);
        }
        if (t%10==0 || t==1) {
            double time = t*dx*dx*0.25;
            printtofile(time,u,n);
        }
    }
}

int main(){
    //Declaring variables
    int n;

    //cout << "Number of gridpoints: ";
    //cin >> n;
    //cout << "Filename to write result too: ";
    //cin >> outfilename;
    double alpha, dx, dt, t_steps;
    n = 10;
    dx = 0.1;
    dt = dx*dx*0.25;
    t_steps = 1000;
    alpha = dt/(dx*dx);
    clock_t start, finish;
    start = clock();

    forward_Euler(n,t_steps,alpha,dx ,dt);

    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;
    finalizePrint(seconds);

    start = clock();
    backwards_Euler(n,t_steps,alpha,dx, dt);
    finish = clock();
    t = ((finish-start));
    seconds = t/CLOCKS_PER_SEC;
    finalizePrint(seconds);


    start = clock();
    crank_Nicolson(n,t_steps,alpha/2,dx, dt);

    finish =clock();
    t = ((finish-start));
    seconds = t/CLOCKS_PER_SEC;

    finalizePrint(seconds);
    string outfilename_ana = "analytic.txt";
    ofile.open(outfilename_ana);
    analytic_Solution(dx,dt,n,t_steps);
    ofile.close();
    return 0;
}
