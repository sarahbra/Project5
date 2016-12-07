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

void GaussElim(int a, vec b, int b_value, int c, int n, vec u, vec v){
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
    for (int i=0; i<=n+1; i++){
        ofile << u(i) << ",";
    }
    ofile << endl;
}

void forward_Euler(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    a_value = c_value = alpha;
    b_value = 1 - 2*alpha;
    double test;
    vec u(n+1);
    vec unew(n+1);
    vec b(n+1);
    vec x = zeros<vec>(n+1);

    u(0) = unew(0) = u(n) = unew(n) = 0.0;
    b(0) = b(0) = b_value;

    //making the vectors
    for (int k=1; k<n; k++){
        b(k) = b_value;
        test = func(dx, k);
        u(k) = test;
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
    vec b(n+1);
    //making the vectors
    for (int k=0; k<=n; k++){
        b(k) = b_value;
    }
}

void crank_Nicholson(int n, int t_steps, double alpha, double dx) {
    double a_value, c_value, b_value;
    vec b(n+1);
    a_value = c_value = -alpha;
    b_value = 2 + 2*alpha;

    //making the vectors
    for (int k=0; k<=n; k++){
        b(k) = b_value;

    //GaussElim(a_value,b, b_value,c_value,u,n,v);

    }
}

int main(){

    //Declaring variables
    int n, a_value, b_value, c_value;

    //cout << "Number of gridpoints: ";
    //cin >> n;
    //cout << "Filename to write result too: ";
    //cin >> outfilename;
    double alpha, dx, dt;
    string outfilename;
    outfilename = "test.txt";
    n = 100;
    dx = 0.1;
    dt = dx*dx*0.25;
    alpha = dt/(dx*dx);
    ofile.open(outfilename);
    clock_t start, finish;
    start = clock();

    forward_Euler(n, 10, alpha, dx);

    finish =clock();
    double t = ((finish-start));
    double seconds = t/CLOCKS_PER_SEC;



    ofile.close();
    return 0;

}
