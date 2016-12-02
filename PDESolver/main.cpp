#include <ctime>
#include <iostream>
#include <cmath>
#include <fstream>
#include "time.h"
#include<string>
using namespace std;

ofstream ofile;

void GaussElim(int a,double b[], int b_value, int c,double d[], int n, double v[]){
    //Forward Substitution
    double m;
    //double v[n+1];// = new double[n+1];
    for (int k=2; k<=n; k++) {
        m = a/b[k-1];
        b[k] = b_value - m*c;
        d[k] = d[k] - m*d[k-1];
    }

    //Backward Substitution
    v[n]= d[n]/b[n];
    cout << "|||||| x[n]" << v[n] << endl;
    for (int k= n-1; k>0; k--) {
        v[k] = (1.0/b[k])*(d[k] - c*v[k+1]);
    }

    v[0] = 0;
    v[n+1] = 0;
}

void function(int n, double ** u_i[], double dx, double i){
    for (int k=0;k<=i+1; k++) {
            u_i(i) = sin(pi*dx*i);
    }
}

void forward_Euler(double alpha) {
    a_value = c_value = alpha;
    b_value = 1 - 2*alpha;
}

void backwards_Euler(double alpha) {
    a_value = c_value = -alpha;
    b_value = 1 + 2*alpha;
}

void crank_Nicholson(double alpha) {
    a_value = c_value = -alpha;
    b_value = 2 + 2*alpha;
}

int main(){

    //Declaring variables
    int n, a_value, b_value, c_value;

    cout << "Number of gridpoints: ";
    cin >> n;
    //cout << "Filename to write result too: ";
    //cin >> outfilename;
    double h, a[n], b[n], c[n], points[n+1];

    //making the vectors
    for (int k=0; k<=n; k++){
        b[k] = b_value;
    }

    for (int k=0; k<=n+1; k++){
        points[k] = k*h;
    }

    double f[n+1];
    double v[n+1];

    function(n, f);

    clock_t start, finish;
    start = clock();

    GaussElim(a_value,b, b_value,c_value, f,n,v);

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
