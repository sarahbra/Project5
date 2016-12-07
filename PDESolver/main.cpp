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

void GaussElim(double a, vec b, double b_value, double c, int n, vec u, vec v){
    //Forward Substitution
    double m;
    for (int k=2; k<=n; k++) {
        m = a/b(k-1);
        b(k) = b_value - m*c;
        u(k) = u(k) - m*u(k-1);
    }

    //Backward Substitution
    v(n)= u(n)/b(n);
    cout << "|||||| x[n]" << v(n) << endl;
    for (int k= n-1; k>0; k--) {
        v(k) = (1.0/b(k))*(u(k) - c*v(k+1));
    }

    v(0) = 0;
    v(n+1) = 0;
}


double func(double dx, double step){
    //double u_i;
    double pi  =3.141592653589793238463;
    return sin(pi*dx*step);

}

void printtofile(int t, vec u, int n ){
    ofile << t << ",";
    for (int i=0; i<=n; i++){
        ofile << u(i) << ",";
    }
    ofile << endl;

}

void forward_Euler(int n, int t_steps, double alpha, double dx) {
    vec u(n+1);
    vec unew(n+1);
    vec b(n+1);

    u(0) = unew(0) = u(n) = unew(n) = 0.0;

    //making the vectors
    for (int k=1; k<n; k++){
        b(k) = b_value;
        u(k) = func(dx, k);
        unew(k) = 0;
    }

    for (int t=1; t<t_steps; t++) {
        for (int i=1; i<n; i++) {
            unew(i) = alpha*u(i-1) + (1-2*alpha) * u(i) + alpha*u(i+1);

        }
        u = unew;
        printtofile(t, u, n);
    }
}

void backwards_Euler(int n, int t_steps, double alpha, double dx) {
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
        v(k) = u(k) = func(dx,k);
    }
    //Implementing boundary conditions
    u(n) = v(n) = u(0) = v(0) = 0.0;

    for (int t=1; t<=t_steps; t++) {
        GaussElim(a_value, b, b_value, c_value, n, v, u);
        u(0) = 0.0;
        u(n) = 0.0;
        for (int i=0; i <= n; i++) {
            v(i) = u(i);
        }
        printtofile(t,u,n);
    }
}

void crank_Nicholson(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    vec b(n+1);
    a_value = c_value = -alpha;
    b_value = 2 + 2*alpha;

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
    for (int t=1; t<=t_steps; t++) {
        for (int i=1; i<n; i++) {
            v(i) = alpha*u(i-1) + (2-2*alpha)*u(i) + alpha*u(i+1);
        }
        v(0) = v(n) = 0;
        GaussElim(a_value,b, b_value,c_value,n+1,v,u);
        u(0) = u(n) = 0.0;
        printtofile(t,u,n);
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
    char* outfilename;
    outfilename = "test.txt";
    n = 10;
    dx = 0.1;
    dt = dx*dx*0.5;
    t_steps = n*n;
    alpha = dt/(dx*dx);
    ofile.open(outfilename);
    clock_t start, finish;
    start = clock();

    //forward_Euler(n, t_steps, alpha, dx);
    backwards_Euler(n,t_steps,alpha,dx);

    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;

    ofile.close();
    return 0;
}
