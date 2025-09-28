#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

//PATRICK CASTILLO 4-775-462    1EE131

// Función para verificar si la matriz es diagonalmente dominante
bool esDiagonalmenteDominante(const vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        double sumaFila = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) sumaFila += fabs(A[i][j]);
        }
        if (fabs(A[i][i]) <= sumaFila) {
            return false; // No es dominante en esta fila
        }
    }
    return true;
}

int main() {
    int n;
    cout << "Ingrese el orden de la matriz (n): ";
    cin >> n;

    // Declaración de la matriz A y el vector b
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    cout << "\nIngrese los coeficientes de la matriz A fila por fila:\n";
    for (int i = 0; i < n; i++) {
        cout << "Fila " << i + 1 << ": ";
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "\nIngrese los valores del vector b:\n";
    for (int i = 0; i < n; i++) {
        cout << "b[" << i + 1 << "]: ";
        cin >> b[i];
    }

    // Aproximación inicial
    vector<double> x(n, 0.0); // Por defecto cero
    cout << "\nDesea introducir una aproximacion inicial? (s/n): ";
    char opcion;
    cin >> opcion;
    if (opcion == 's' || opcion == 'S') {
        for (int i = 0; i < n; i++) {
            cout << "x0[" << i + 1 << "]: ";
            cin >> x[i];
        }
    }

    double tol;
    int maxIter;
    cout << "\nIngrese la tolerancia: ";
    cin >> tol;
    cout << "Ingrese el numero maximo de iteraciones: ";
    cin >> maxIter;

    // Verificar si la matriz es diagonalmente dominante
    if (!esDiagonalmenteDominante(A)) {
        cout << "\nADVERTENCIA: La matriz no es diagonalmente dominante.\n";
        cout << "El metodo de Jacobi puede NO converger.\n";
    }

    // Variables para iterar
    vector<double> xNuevo(n, 0.0);
    vector<double> error(n, 100.0);

    cout << fixed << setprecision(6);
    cout << "\nTabla de iteraciones (Jacobi):\n";
    cout << "Iter\t";
    for (int i = 0; i < n; i++) {
        cout << "X" << i + 1 << "\tError%\t";
    }
    cout << endl;

    bool convergencia = false;
//PATRICK CASTILLO 4-775-462    1EE131
    for (int iter = 1; iter <= maxIter; iter++) {
        // Calcular nueva aproximación
        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (i != j) suma -= A[i][j] * x[j];
            }
            xNuevo[i] = suma / A[i][i];
        }

        // Calcular errores relativos
        convergencia = true;
        for (int i = 0; i < n; i++) {
            if (fabs(xNuevo[i]) > 1e-12) { // evitar división por cero
                error[i] = fabs((xNuevo[i] - x[i]) / xNuevo[i]) * 100.0;
            } else {
                error[i] = fabs(xNuevo[i] - x[i]);
            }
            if (error[i] > tol * 100.0) { // convertir tol a porcentaje
                convergencia = false;
            }
        }

        // Imprimir fila de la tabla
        cout << iter << "\t";
        for (int i = 0; i < n; i++) {
            cout << xNuevo[i] << "\t" << error[i] << "\t";
        }
        cout << endl;
//PATRICK CASTILLO 4-775-462    1EE131
        // Actualizar x
        x = xNuevo;

        if (convergencia) {
            cout << "\nEl metodo ha convergido en la iteracion " << iter << ".\n";
            break;
        }
    }

    if (!convergencia) {
        cout << "\nEl metodo no convergio en " << maxIter << " iteraciones.\n";
    }

    return 0;
}
