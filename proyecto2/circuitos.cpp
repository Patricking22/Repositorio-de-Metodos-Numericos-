#include <iostream> 
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

/* ================================================================
   Programa: Resolución de sistemas de ecuaciones lineales
   Aplicación: Análisis de mallas o nodos eléctricos
   Autor: Estudiante de Ingeniería Eléctrica
   Lenguaje: C++
   ================================================================
*/

// ----------- Declaración de funciones ----------------
void metodoJacobi(vector<vector<double>> A, vector<double> b, double tol, int maxIter);
void metodoGaussSeidel(vector<vector<double>> A, vector<double> b, double tol, int maxIter);
void eliminacionGauss(vector<vector<double>> A, vector<double> b);
void doolittleLU(vector<vector<double>> A, vector<double> b);
void matrizInversa(vector<vector<double>> A, vector<double> b);

// ================================================================
// FUNCION PRINCIPAL
// ================================================================
int main() {
    int n;
    char repetir;

    do {
        cout <<"\n-------UNIVERSIDAD TECNOLOGICA DE PANAMA------\n";
        cout <<"\n---FACULTAD DE INGENIERIA ELECTRICA--- \n";
        cout << "\n--INGENIERIA ELECTRICA Y ELECTRONICA--\n";    
        cout <<"\n -PROYECTO 2 - METODOS NUMERICOS- \n";   

        cout << "\n=====================================================\n";
        cout << "   RESOLUCION DE SISTEMAS DE ECUACIONES LINEALES\n";
        cout << "   Aplicacion: Analisis de Mallas o Nodos Electricos\n";
        cout << "=====================================================\n";

        // --- Nueva sección: lista descriptiva de métodos ---
        cout << "\nLISTA DE METODOS DISPONIBLES:\n";
        cout << "-----------------------------------------------------\n";
        cout << "1. Metodo de Jacobi (Iterativo):\n";
        cout << "   Calcula cada variable usando los valores anteriores,\n";
        cout << "   adecuado para matrices diagonales dominantes.\n\n";
        cout << "2. Metodo de Gauss-Seidel (Iterativo):\n";
        cout << "   Similar a Jacobi, pero actualiza cada variable\n";
        cout << "   inmediatamente con el nuevo valor, lo que acelera la convergencia.\n\n";
        cout << "3. Eliminacion de Gauss (Directo):\n";
        cout << "   Reduce el sistema a una forma triangular y aplica sustitucion regresiva.\n\n";
        cout << "4. Metodo de Doolittle (LU):\n";
        cout << "   Descompone la matriz A en L (inferior) y U (superior),\n";
        cout << "   ideal para resolver varios sistemas con la misma A.\n\n";
        cout << "5. Metodo de la Matriz Inversa (Directo):\n";
        cout << "   Calcula la matriz inversa de A y obtiene x = A^-1 * b.\n";
        cout << "-----------------------------------------------------\n";

        cout << "\nIngrese el orden del sistema (n): ";
        cin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);

        cout << "\nIngrese los coeficientes de la matriz A (" << n << "x" << n << "):\n";
        for (int i = 0; i < n; i++) {
            cout << "Fila " << i + 1 << ": ";
            for (int j = 0; j < n; j++) cin >> A[i][j];
        }

        cout << "\nIngrese el vector b (voltajes o terminos independientes):\n";
        for (int i = 0; i < n; i++) {
            cout << "b[" << i + 1 << "]: ";
            cin >> b[i];
        }

        cout << "\nSeleccione el metodo de resolucion:\n";
        cout << "1. Metodo de Jacobi (Iterativo)\n";
        cout << "2. Metodo de Gauss-Seidel (Iterativo)\n";
        cout << "3. Eliminacion de Gauss (Directo)\n";
        cout << "4. Metodo de Doolittle (LU)\n";
        cout << "5. Metodo de la Matriz Inversa (Directo)\n";
        cout << "Opcion: ";
        int opcion;
        cin >> opcion;

        double tol = 0.0001;
        int maxIter = 20;

        if (opcion == 1 || opcion == 2) {
            cout << "\nIngrese la tolerancia: ";
            cin >> tol;
            cout << "Ingrese el numero maximo de iteraciones: ";
            cin >> maxIter;
        }

        cout << "\n-----------------------------------------------------\n";

        switch (opcion) {
            case 1: metodoJacobi(A, b, tol, maxIter); break;
            case 2: metodoGaussSeidel(A, b, tol, maxIter); break;
            case 3: eliminacionGauss(A, b); break;
            case 4: doolittleLU(A, b); break;
            case 5: matrizInversa(A, b); break;
            default: cout << "Opcion no valida.\n"; break;
        }

        cout << "-----------------------------------------------------\n";
        cout << "\nDesea resolver otro sistema? (S/N): ";
        cin >> repetir;
    } while (repetir == 'S' || repetir == 's');

    cout << "\nPrograma finalizado. Gracias por usar el sistema.\n";
    return 0;
}

// ================================================================
// METODO DE JACOBI
// ================================================================
void metodoJacobi(vector<vector<double>> A, vector<double> b, double tol, int maxIter) {
    int n = A.size();
    vector<double> x(n, 0.0), xNuevo(n, 0.0);

    cout << "\n*** Metodo de Jacobi ***\n";
    cout << left << setw(8) << "Iter";
    for (int i = 0; i < n; i++) cout << setw(12) << ("x" + to_string(i+1));
    cout << setw(12) << "Error (%)" << endl;

    for (int iter = 1; iter <= maxIter; iter++) {
        for (int i = 0; i < n; i++) {
            double suma = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) suma += A[i][j] * x[j];
            }
            xNuevo[i] = (b[i] - suma) / A[i][i];
        }

        double errorMax = 0.0;
        for (int i = 0; i < n; i++) {
            double error = fabs((xNuevo[i] - x[i]) / (xNuevo[i] + 1e-12)) * 100.0;
            if (error > errorMax) errorMax = error;
        }

        cout << setw(8) << iter;
        for (int i = 0; i < n; i++) cout << setw(12) << fixed << setprecision(6) << xNuevo[i];
        cout << setw(12) << errorMax << endl;

        if (errorMax < tol * 100) {
            cout << "\nConvergencia alcanzada en " << iter << " iteraciones.\n";
            break;
        }
        x = xNuevo;
    }

    cout << "\nSolucion aproximada:" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << x[i] << endl;
}

// ================================================================
// METODO DE GAUSS-SEIDEL
// ================================================================
void metodoGaussSeidel(vector<vector<double>> A, vector<double> b, double tol, int maxIter) {
    int n = A.size();
    vector<double> x(n, 0.0);

    cout << "\n*** Metodo de Gauss-Seidel ***\n";
    cout << left << setw(8) << "Iter";
    for (int i = 0; i < n; i++) cout << setw(12) << ("x" + to_string(i+1));
    cout << setw(12) << "Error (%)" << endl;

    for (int iter = 1; iter <= maxIter; iter++) {
        vector<double> xViejo = x;
        for (int i = 0; i < n; i++) {
            double suma = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) suma += A[i][j] * x[j];
            }
            x[i] = (b[i] - suma) / A[i][i];
        }

        double errorMax = 0.0;
        for (int i = 0; i < n; i++) {
            double error = fabs((x[i] - xViejo[i]) / (x[i] + 1e-12)) * 100.0;
            if (error > errorMax) errorMax = error;
        }

        cout << setw(8) << iter;
        for (int i = 0; i < n; i++) cout << setw(12) << fixed << setprecision(6) << x[i];
        cout << setw(12) << errorMax << endl;

        if (errorMax < tol * 100) {
            cout << "\nConvergencia alcanzada en " << iter << " iteraciones.\n";
            break;
        }
    }

    cout << "\nSolucion aproximada:" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << x[i] << endl;
}

// ================================================================
// ELIMINACION DE GAUSS
// ================================================================
void eliminacionGauss(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<double> x(n);

    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        double suma = 0.0;
        for (int j = i + 1; j < n; j++) suma += A[i][j] * x[j];
        x[i] = (b[i] - suma) / A[i][i];
    }

    cout << "\n*** Eliminacion de Gauss ***" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
}

// ================================================================
// DOOLITTLE (LU)
// ================================================================
void doolittleLU(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0.0)), U(n, vector<double>(n, 0.0));
    vector<double> y(n), x(n);

    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
        for (int j = i; j < n; j++) {
            double suma = 0.0;
            for (int k = 0; k < i; k++) suma += L[i][k] * U[k][j];
            U[i][j] = A[i][j] - suma;
        }
        for (int j = i + 1; j < n; j++) {
            double suma = 0.0;
            for (int k = 0; k < i; k++) suma += L[j][k] * U[k][i];
            L[j][i] = (A[j][i] - suma) / U[i][i];
        }
    }

    y[0] = b[0];
    for (int i = 1; i < n; i++) {
        double suma = 0.0;
        for (int j = 0; j < i; j++) suma += L[i][j] * y[j];
        y[i] = b[i] - suma;
    }

    for (int i = n - 1; i >= 0; i--) {
        double suma = 0.0;
        for (int j = i + 1; j < n; j++) suma += U[i][j] * x[j];
        x[i] = (y[i] - suma) / U[i][i];
    }

    cout << "\n*** Metodo de Doolittle (LU) ***" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
}

// ================================================================
// MATRIZ INVERSA
// ================================================================
void matrizInversa(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> inv(n, vector<double>(n, 0.0));

    double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    inv[0][0] = A[1][1] / det;
    inv[0][1] = -A[0][1] / det;
    inv[1][0] = -A[1][0] / det;
    inv[1][1] = A[0][0] / det;

    vector<double> x(n);
    for (int i = 0; i < n; i++) {
        x[i] = inv[i][0]*b[0] + inv[i][1]*b[1];
    }

    cout << "\n*** Metodo de la Matriz Inversa ***" << endl;
    for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << fixed << setprecision(6) << x[i] << endl;
}
