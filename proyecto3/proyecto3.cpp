#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace std;

double f_integrand(double x) {
    // f(x) = (3x^3 - 24x^2 + 48x + 5) / (x^2 - 8x + 16)
    double num = 3.0*x*x*x - 24.0*x*x + 48.0*x + 5.0;
    double den = x*x - 8.0*x + 16.0; // (x-4)^2, no cero en [-3,2]
    return num / den;
}

// ---------- Métodos de integración ----------
void metodoTrapecio(double a, double b, int n) {
    cout << "\n--- Metodo del Trapecio Compuesto ---\n";
    cout << "Intervalo: [" << a << ", " << b << "], n = " << n << "\n";
    double h = (b - a) / n;
    vector<double> xi(n+1), yi(n+1);
    for (int i = 0; i <= n; ++i) {
        xi[i] = a + i*h;
        yi[i] = f_integrand(xi[i]);
    }

    cout << fixed << setprecision(6);
    cout << "\ni\t  xi\t\t  f(xi)\n";
    for (int i = 0; i <= n; ++i) {
        cout << setw(2) << i << "   " << setw(10) << xi[i] << "   " << setw(10) << yi[i] << "\n";
    }

    // Suma para trapecio
    double suma = 0.0;
    for (int i = 1; i < n; ++i) suma += yi[i];
    double integral = (h/2.0) * (yi[0] + 2.0*suma + yi[n]);
    cout << "\nSumas intermedias:\n";
    cout << "h = " << h << "\n";
    cout << "Sum interiores (sum_{i=1}^{n-1} f(xi)) = " << suma << "\n";
    cout << "Integral (Trapecio) = (h/2) * (f0 + 2*sum + fn) = " << integral << "\n";
}

void metodoSimpson13(double a, double b, int n) {
    cout << "\n--- Metodo de Simpson 1/3 Compuesto ---\n";
    cout << "Intervalo: [" << a << ", " << b << "], n = " << n << " (n debe ser par)\n";
    if (n % 2 == 1) {
        cout << "Error: n debe ser par para Simpson 1/3. Ajustando n = n+1.\n";
        ++n;
    }
    double h = (b - a) / n;
    vector<double> xi(n+1), yi(n+1);
    for (int i = 0; i <= n; ++i) {
        xi[i] = a + i*h;
        yi[i] = f_integrand(xi[i]);
    }

    cout << fixed << setprecision(6);
    cout << "\ni\t  xi\t\t  f(xi)\n";
    for (int i = 0; i <= n; ++i) cout << setw(2) << i << "   " << setw(10) << xi[i] << "   " << setw(10) << yi[i] << "\n";

    double suma_odd = 0.0, suma_even = 0.0;
    for (int i = 1; i < n; ++i) {
        if (i % 2 == 1) suma_odd += yi[i];
        else suma_even += yi[i];
    }
    double integral = (h/3.0) * (yi[0] + 4.0*suma_odd + 2.0*suma_even + yi[n]);
    cout << "\nSumas intermedias:\n";
    cout << "h = " << h << "\n";
    cout << "Sum odd (i impar) = " << suma_odd << "\n";
    cout << "Sum even (i par) = " << suma_even << "\n";
    cout << "Integral (Simpson 1/3) = " << integral << "\n";
}

void metodoSimpson38(double a, double b, int n) {
    cout << "\n--- Metodo de Simpson 3/8 Compuesto ---\n";
    cout << "Intervalo: [" << a << ", " << b << "], n = " << n << " (n debe ser múltiplo de 3)\n";
    if (n % 3 != 0) {
        int n_old = n;
        n = n + (3 - (n % 3));
        cout << "Ajustando n de " << n_old << " a " << n << " (múltiplo de 3) para Simpson 3/8.\n";
    }
    double h = (b - a) / n;
    vector<double> xi(n+1), yi(n+1);
    for (int i = 0; i <= n; ++i) {
        xi[i] = a + i*h;
        yi[i] = f_integrand(xi[i]);
    }

    cout << fixed << setprecision(6);
    cout << "\ni\t  xi\t\t  f(xi)\n";
    for (int i = 0; i <= n; ++i) cout << setw(2) << i << "   " << setw(10) << xi[i] << "   " << setw(10) << yi[i] << "\n";

    double sum1 = 0.0; // indices not multiple of 3
    double sum3 = 0.0; 
    for (int i = 1; i < n; ++i) {
        if (i % 3 == 0) sum3 += yi[i];
        else sum1 += yi[i];
    }
    double integral = (3.0*h/8.0) * (yi[0] + 3.0*sum1 + 2.0*sum3 + yi[n]);
    cout << "\nSumas intermedias:\n";
    cout << "h = " << h << "\n";
    cout << "Sum (i%3 != 0) = " << sum1 << "\n";
    cout << "Sum (i%3 == 0, i interior) = " << sum3 << "\n";
    cout << "Integral (Simpson 3/8) = " << integral << "\n";
}


double lagrangeEvaluate(const vector<double>& xs, const vector<double>& ys, double x0) {
    int n = xs.size();
    double px = 0.0;
    cout << fixed << setprecision(8);
    for (int i = 0; i < n; ++i) {
        double Li = 1.0;
        cout << "\nTermino L_" << i << "(x): (";
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            cout << "(x - " << xs[j] << ")/(" << xs[i] << " - " << xs[j] << ")";
            Li *= (x0 - xs[j]) / (xs[i] - xs[j]);
            if (j != n-1 && !(j==n-2 && i==n-1)) cout << " * ";
        }
        double contrib = ys[i] * Li;
        cout << ") evaluado en x0 = " << x0 << " da L_" << i << "(" << x0 << ") = " << Li
             << " contribucion = " << contrib << "\n";
        px += contrib;
    }
    return px;
}


double newtonDividedDifferences(const vector<double>& xs, const vector<double>& ys, double x0) {
    int n = xs.size();
    vector<vector<double>> dd(n, vector<double>(n, 0.0));
    // Primera columna = ys
    for (int i = 0; i < n; ++i) dd[i][0] = ys[i];

    // Construir tabla
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            double num = dd[i+1][j-1] - dd[i][j-1];
            double den = xs[i+j] - xs[i];
            dd[i][j] = num / den;
        }
    }

    // Mostrar tabla
    cout << fixed << setprecision(10);
    cout << "\nTabla de diferencias divididas (formato dd[i][j] donde j es orden):\n";
    for (int i = 0; i < n; ++i) {
        cout << "i=" << i << ": ";
        for (int j = 0; j < n - i; ++j) cout << setw(14) << dd[i][j];
        cout << "\n";
    }

    // Evaluar polinomio de Newton en x0
    double px = dd[0][0];
    cout << "\nEvaluacion polinomio de Newton:\n";
    cout << "Term 0: " << dd[0][0] << "\n";
    double prod = 1.0;
    for (int k = 1; k < n; ++k) {
        prod *= (x0 - xs[k-1]);
        double term = dd[0][k] * prod;
        cout << "Term " << k << ": dd[0]["<<k<<"] = " << dd[0][k] << "  prod = " << prod << "  contrib = " << term << "\n";
        px += term;
    }
    return px;
}

// Extrapolacion polinomial grado 3 usando últimos 4 puntos (xs[n-4..n-1])
double extrapolacionPolinomialDeg3(const vector<double>& xs_full, const vector<double>& ys_full, double x0) {
    int nfull = xs_full.size();
    if (nfull < 4) {
        cerr << "No hay suficientes puntos para extrapolacion grado 3.\n";
        return NAN;
    }
    vector<double> xs(4), ys(4);
    for (int i = 0; i < 4; ++i) {
        xs[i] = xs_full[nfull - 4 + i];
        ys[i] = ys_full[nfull - 4 + i];
    }
    cout << "\nUsando puntos para extrapolacion (grado 3):\n";
    for (int i = 0; i < 4; ++i) cout << "x=" << xs[i] << "  f=" << ys[i] << "\n";

    // Usamos Newton  sobre esos 4 puntos (grado 3)
    double px = newtonDividedDifferences(xs, ys, x0);
    return px;
}

int main() {
    cout << setprecision(10);
    cout << "Programa: Integracion numerica e interpolacion/extrapolacion\n";
    cout << "Opciones:\n";
    cout << "1) Integracion numerica (Trapecio, Simpson 1/3, Simpson 3/8)\n";
    cout << "2) Interpolacion / Extrapolacion (Extrapolacion polinomial, Lagrange, Newton)\n";
    cout << "Elija 1 o 2: ";
    int modo; cin >> modo;

    if (modo == 1) {
        // Integracion del problema dado
        double a = -3.0, b = 2.0;
        cout << "\nHas elegido INTEGRACION numerica para:\n";
        cout << "Integral desde " << a << " hasta " << b << " de (3x^3 -24x^2 +48x +5)/(x^2 -8x +16) dx\n";

        int n;
        cout << "\nMETODO DEL TRAPECIO: ingrese numero de subintervalos n (>=1): ";
        cin >> n;
        if (n < 1) n = 1;
        metodoTrapecio(a, b, n);

        cout << "\nMETODO DE SIMPSON 1/3: ingrese numero de subintervalos n (par): ";
        cin >> n;
        if (n < 2) n = 2;
        metodoSimpson13(a, b, n);

        cout << "\nMETODO DE SIMPSON 3/8: ingrese numero de subintervalos n (múltiplo de 3): ";
        cin >> n;
        if (n < 3) n = 3;
        metodoSimpson38(a, b, n);

        cout << "\nFin de calculos de integracion.\n";
    }
    else if (modo == 2) {
        // Interpolacion / extrapolacion: tabla dada
        vector<double> xs = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
        vector<double> ys = {0.3, 0.31, 0.32, 0.33, 0.34, 0.45, 0.46, 0.47};
        cout << "\nHas elegido INTERPOLACION / EXTRAPOLACION.\n";
        cout << "Tabla de datos (mostrar):\n";
        cout << "i\t x\t f(x)\n";
        for (size_t i = 0; i < xs.size(); ++i) cout << i << "\t" << xs[i] << "\t" << ys[i] << "\n";

        double x0 = 3.5;
        cout << "\nQueremos estimar f(" << x0 << ").\n";

        // Metodo 4: Extrapolacion polinomial (grado 3 usando ultimos 4 puntos)
        cout << "\n================== Metodo 4: Extrapolacion Polinomial (grado 3) ==================\n";
        double p_extrap = extrapolacionPolinomialDeg3(xs, ys, x0);
        cout << "\nResultado (Extrapolacion polinomial deg 3): f(" << x0 << ") = " << p_extrap << "\n";

        // Metodo 5: Lagrange (usar todos los puntos)
        cout << "\n================== Metodo 5: Lagrange (8 puntos) ==================\n";
        double p_lagr = lagrangeEvaluate(xs, ys, x0);
        cout << "\nResultado (Lagrange con 8 puntos): f(" << x0 << ") = " << p_lagr << "\n";

        // Metodo 6: Newton (dif. divididas) con todos los puntos
        cout << "\n================== Metodo 6: Newton (Diferencias Divididas, 8 puntos) ==================\n";
        double p_newton = newtonDividedDifferences(xs, ys, x0);
        cout << "\nResultado (Newton con 8 puntos): f(" << x0 << ") = " << p_newton << "\n";

        cout << "\nFin de calculos de interpolacion/extrapolacion.\n";
    }
    else {
        cout << "Opcion invalida. Finalizando.\n";
    }

    cout << "\nPrograma terminado.\n";
    return 0;
}
