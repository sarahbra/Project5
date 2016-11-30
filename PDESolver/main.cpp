#include <ctime>
#include <iostream>
using namespace std;

void fillVec(double[], int, string);

int main()
{
    //Defining variables
    int n;
    cout << "Dimension of matrix:";
    cin >> n;
    double a[n-1], b[n], c[n-1], vecu[n], vecf[n];

    //Assigning random numbers to vector elements (other than zero, making the tridiagonal matrix unsolvable)
    srand((unsigned)time(NULL));
    fillVec(b,n,"b");
    fillVec(a,n-1,"a");
    fillVec(c,n-1,"c");
    fillVec(vecf,n,"f");

    //Forwards substitution
    for (int k=1; k<=n; k++) {
        b[k] = b[k] - (a[k-1]*c[k-1])/b[k-1];
        vecf[k] = vecf[k] - (a[k-1]*vecf[k-1])/b[k-1];
    }

    //Backwards substitution
    vecu[n] = b[n];
    double temp;
    cout << "u" << endl;
    for (int k=n-1; k>=0; k--) {
        temp = c[k]*vecu[k+1];
        vecu[k] = 1.0/b[k]*(vecf[k] - temp);
        cout << vecu[k] << endl;
    }
    return 0;
}

void fillVec(double arr[], int n, string nb)  {
    cout << nb << endl;
    for (int i=0; i<n; i++) {
        arr[i] = rand()%100 + 1;
        cout << arr[i] << endl;
    }
}
