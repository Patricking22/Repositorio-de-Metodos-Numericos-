#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

//PATRICK CASTILLO 4-775-462    1EE131
double evalPoli(const vector<double>& coef, double x) {
double result = 0.0;
int grado = coef.size() - 1;
for (int i = 0; i <= grado; i++) {
result += coef[i] * pow(x, grado - i);
}
return result;
}


vector<double> derivarPolinomio(const vector<double>& coef) {
vector<double> derivada;
int grado = coef.size() - 1;
for (int i = 0; i < grado; i++) {
derivada.push_back(coef[i] * (grado - i));
}
return derivada;
}

int main() {
int grado;
cout << "Metodo de Newton y Raphson\n";
cout << "Ingrese el grado del polinomio: ";
cin >> grado;

vector<double> coef(grado + 1);
cout << "Ingrese los coeficientes de mayor a menor grado:\n";
for (int i = 0; i <= grado; i++) {
    cout << "Coeficiente x^" << (grado - i) << ": ";
    cin >> coef[i];
}
//PATRICK CASTILLO 4-775-462    1EE131
// Derivada del polinomio
vector<double> derivada = derivarPolinomio(coef);

double x0;
cout << "Ingrese el valor inicial x0: ";
cin >> x0;

int maxIter;
cout << "Ingrese el numero maximo de iteraciones: ";
cin >> maxIter;

double tolerancia;
cout << "Ingrese la tolerancia de error: ";
cin >> tolerancia;

cout << fixed << setprecision(6);
cout << "\nTabla de iteraciones (Newton-Raphson):\n";
cout << "n\t x_n\t\t f(x_n)\t\t f'(x_n)\t x_{n+1}\t\t Error\n";
cout << "-----------------------------------------------------------------------------\n";

double x = x0;
for (int n = 1; n <= maxIter; n++) {
    double fx = evalPoli(coef, x);
    double fpx = evalPoli(derivada, x);
//PATRICK CASTILLO 4-775-462    1EE131
    if (fabs(fpx) < 1e-12) {
        cout << "La derivada es cero, no se puede continuar.\n";
        break;
    }

    double x1 = x - fx / fpx;
    double error = fabs(x1 - x);

    cout << n << "\t " << x << "\t " << fx << "\t " << fpx << "\t " 
         << x1 << "\t " << error << "\n";

    if (error < tolerancia) {
        cout << "\nSe encontro una raiz aproximada en x = " << x1 << "\n";
        break;
    }
    x = x1;
}

return 0;

}
