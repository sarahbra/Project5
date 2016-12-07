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

void GaussElim(double a, vec b(), double b_value, double c, int n, vec u(), vec v()){
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

double function(double dx, double step){
    u_i = sin(pi*dx*step);
    return u_i;
}

void forward_Euler(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    a_value = c_value = alpha;
    b_value = 1 - 2*alpha;
    vec b(n);
    vec u(n+1);
    vec unew(n+1);

    u(0) = unew(0) = u(n) = unew(n) = 0.0;
    b(0) = b(n) = b_value;

    //making the vectors
    for (int k=1; k<n; k++){
        b(k) = b_value;
        u(k) = function(dx,i);
        unew(k) = 0;
    }

    for (int t=1; t<t_steps; t++) {
        for (int i=1; i<n; i++) {
            unew(i) = alpha*u(i-1) + (1-2*alpha) * u(i) + alpha*u(i+1);
        }
    }
}

void backwards_Euler(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    a_value = c_value = -alpha;
    b_value = 1 + 2*alpha;
    vec b(n);

    vec u(n+1);
    vec v(n+1);

    //making the vectors
    for (int k=0; k<=n; k++){
        b(k) = b_value;
    }
    for (int k=1; k<n; k++) {
        v(i) = u(i) = function(dx,i);
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
    }
}

void crank_Nicholson(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    a_value = c_value = -alpha;
    b_value = 2 + 2*alpha;
    vec b(n);

    vec u(n+1);
    vec v(n+1);

    //making the vectors
    for (int k=1; k<n; k++){
        b(k) = b_value;
        u(k) = function(dx,k);
    }

    b(0) = b(n) = b_value;
    u(0) = u(n) = 0.0;

    for (int t=1; t<=t_steps; t++) {
        for (int i=1; i<n; i++) {
            v(i) = alpha*u(i-1) + (2-2*alpha)*u(i) + alpha*u(i+1);
        }
        v(0) = v(n) = 0;
        GaussElim(a_value,b, b_value,c_value,n+1,v,u);
        u(0) = u(n) = 0.0;
    }
}

int main(){
    cout << "Number of gridpoints: ";
    cin >> n;
    //cout << "Filename to write result too: ";
    //cin >> outfilename;
    double alpha, dx, dt;

    dx = 1.0/10;
    dt = dx*dx*0.25;
    t_steps = n^2;
    alpha = dt/(dx*dx);

    clock_t start, finish;
    start = clock();

    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;
    string outfilename;

    cout << "Please enter a file name to write: ";
    cin >> outfilename;

    ofile.open(outfilename);
    int i;
    ofile << "x:  " <<"computed_derivative" << "n = "<< n << "runtime: "<< seconds << endl;
    for (i=0; i<=n+1; i++){
        ofile << points[i] << "   " << v[i] << "\n";
    }
    ofile.close();
    return 0;
}
