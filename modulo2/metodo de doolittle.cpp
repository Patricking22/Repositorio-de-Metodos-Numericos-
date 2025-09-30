#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

int main() {
int n;
cout << "Metodo de Doolittle (Descomposicion LU y resolucion de Ax=b)\n";
cout << "Ingrese el orden de la matriz: ";
cin >> n;


vector<vector<double>> A(n, vector<double>(n));
vector<vector<double>> L(n, vector<double>(n, 0));
vector<vector<double>> U(n, vector<double>(n, 0));
vector<double> b(n), y(n), x(n);

cout << "Ingrese los elementos de la matriz A (" << n << "x" << n << "):\n";
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        cout << "A[" << i+1 << "][" << j+1 << "] = ";
        cin >> A[i][j];
    }
}

cout << "Ingrese el vector b (" << n << " elementos):\n";
for (int i = 0; i < n; i++) {
    cout << "b[" << i+1 << "] = ";
    cin >> b[i];
}

// Doolittle: diagonal de L = 1
for (int i = 0; i < n; i++) {
    L[i][i] = 1;
}

// Calculo de L y U
for (int k = 0; k < n; k++) {
    // Fila de U
    for (int j = k; j < n; j++) {
        double suma = 0;
        for (int p = 0; p < k; p++) {
            suma += L[k][p] * U[p][j];
        }
        U[k][j] = A[k][j] - suma;
    }

    // Columna de L
    for (int i = k+1; i < n; i++) {
        double suma = 0;
        for (int p = 0; p < k; p++) {
            suma += L[i][p] * U[p][k];
        }
        if (U[k][k] == 0) {
            cout << "Error: pivote nulo.\n";
            return 1;
        }
        L[i][k] = (A[i][k] - suma) / U[k][k];
    }
}

cout << fixed << setprecision(3);

cout << "\nMatriz L:\n";
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        cout << setw(8) << L[i][j] << " ";
    }
    cout << "\n";
}

cout << "\nMatriz U:\n";
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        cout << setw(8) << U[i][j] << " ";
    }
    cout << "\n";
}

// Sustitucion hacia adelante: L * y = b
for (int i = 0; i < n; i++) {
    double suma = 0;
    for (int j = 0; j < i; j++) {
        suma += L[i][j] * y[j];
    }
    y[i] = b[i] - suma;
}

// Sustitucion hacia atras: U * x = y
for (int i = n-1; i >= 0; i--) {
    double suma = 0;
    for (int j = i+1; j < n; j++) {
        suma += U[i][j] * x[j];
    }
    if (U[i][i] == 0) {
        cout << "Error: division por cero en U.\n";
        return 1;
    }
    x[i] = (y[i] - suma) / U[i][i];
}

cout << "\nSolucion del sistema Ax = b:\n";
for (int i = 0; i < n; i++) {
    cout << "x[" << i+1 << "] = " << x[i] << "\n";
}

return 0;


}
