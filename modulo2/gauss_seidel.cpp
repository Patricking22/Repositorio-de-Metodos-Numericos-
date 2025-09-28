#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

//PATRICK  CASTILLO 4-775-462    1EE131
int main() {
    int n;
    cout << "Ingrese el orden de la matriz (n): ";
    cin >> n;

    // Matriz de coeficientes A y vector b
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n);

    cout << "\nIngrese la matriz A (coeficientes):\n";
    for (int i = 0; i < n; i++) {
        cout << "Fila " << i + 1 << ": ";
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "\nIngrese el vector b (terminos independientes):\n";
    for (int i = 0; i < n; i++) {
        cout << "b[" << i + 1 << "]: ";
        cin >> b[i];
    }

    // Aproximación inicial
    vector<double> x(n, 0.0);
    cout << "\nDesea introducir una aproximacion inicial? (s/n): ";
    char opcion;
    cin >> opcion;
    if (opcion == 's' || opcion == 'S') {
        for (int i = 0; i < n; i++) {
            cout << "x0[" << i + 1 << "]: ";
            cin >> x[i];
        }
    }
    //PATRICK  CASTILLO 4-775-462    1EE131
    double tol;
    int maxIter;
    cout << "\nIngrese la tolerancia: ";
    cin >> tol;
    cout << "Ingrese el numero maximo de iteraciones: ";
    cin >> maxIter;

    cout << fixed << setprecision(6);
    cout << "\nTabla de iteraciones (Gauss-Seidel):\n";
    cout << "Iter\t";
    for (int i = 0; i < n; i++) {
        cout << "X" << i + 1 << "\tError%\t";
    }
    cout << endl;

    vector<double> xViejo(n, 0.0);
    vector<double> error(n, 100.0);
    bool convergencia = false;

    for (int iter = 1; iter <= maxIter; iter++) {
        xViejo = x; // guardar valores anteriores

        // Calcular nueva aproximación usando Gauss-Seidel
        for (int i = 0; i < n; i++) {
            double suma = b[i];
            for (int j = 0; j < n; j++) {
                if (i != j) suma -= A[i][j] * x[j]; // aquí x[j] ya puede ser actualizado
            }
            x[i] = suma / A[i][i];
        }
        //PATRICK  CASTILLO 4-775-462    1EE131
        // Calcular errores relativos
        convergencia = true;
        for (int i = 0; i < n; i++) {
            if (fabs(x[i]) > 1e-12)
                error[i] = fabs((x[i] - xViejo[i]) / x[i]) * 100.0;
            else
                error[i] = fabs(x[i] - xViejo[i]);

            if (error[i] > tol * 100.0) convergencia = false;
        }

        // Imprimir resultados de esta iteración
        cout << iter << "\t";
        for (int i = 0; i < n; i++) {
            cout << x[i] << "\t" << error[i] << "\t";
        }
        cout << endl;

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

