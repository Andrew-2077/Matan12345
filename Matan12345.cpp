#include <omp.h>
#include <iostream>
#include <math.h>
#include <Windows.h>

double func(double* y, double time, int i)
{
    double w;
    switch (i) {
    case 0:
        w = y[1];
        break;
    case 1:
        w = -0.01 * y[0] * exp(-time);
        break;
    default:
        break;
    }
    return w;
}

void Euler1()
{
    const int n = 2;
    SetConsoleOutputCP(1251);
    double t0 = 0.0, tmax = 7.0, tau = 0.01, t = t0;
    double y[n] = { 0.0, 0.5 }, yy[n] = { 0.0 }, tn, tk, deltat;
    double a[n][n], b[n], eps;
    tn = omp_get_wtime();
    for (double t = t0; t <= tmax; t += tau)
    {
        for (int i = 0; i < n; i++)
            yy[i] = y[i] + tau * func(y, t, i);
        for (int i = 0; i < n; i++)
            y[i] = yy[i];
    }
    tk = omp_get_wtime();
    deltat = tk - tn;
    printf("Массив параметров: ");
    for (int i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\nВремя вычислений: %f\n\n", deltat);
}

void RK2()
{
    const int n = 2;
    SetConsoleOutputCP(1251);
    double t0 = 0.0, tmax = 7.0, tau = 0.01, t = t0;
    double y[n] = { 0.0, 0.5 }, yy[n] = { 0.0 }, tn, tk, deltat, ff[n] = { 0.0 };
    tn = omp_get_wtime();
    for (double t = t0; t <= tmax; t += tau)
    {
        for (int i = 0; i < n; i++)
            yy[i] = y[i] + (tau * func(y, t, i)) / 2.0;
        for (int i = 0; i < n; i++)
            ff[i] = func(yy, t + 0.5 * tau, i);
        for (int i = 0; i < n; i++)
            y[i] += tau * ff[i];
    }
    tk = omp_get_wtime();
    deltat = tk - tn;
    printf("Массив параметров: ");
    for (int i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\nВремя вычислений: %f\n\n", deltat);
}

void PredKor()
{
    const int n = 2;
    SetConsoleOutputCP(1251);
    double t0 = 0.0, tmax = 7.0, tau = 0.01, t = t0;
    double y[n] = { 0.0, 0.5 }, yy[n] = { 0.0 }, tn, tk, deltat, ff[n] = { 0.0 };
    double a[n] = { 0.0 };
    tn = omp_get_wtime();
    for (double t = t0; t <= tmax; t += tau)
    {
        for (int i = 0; i < n; i++)
            a[i] = func(y, t, i);
        for (int i = 0; i < n; i++)
            yy[i] = y[i] + tau * a[i];
        for (int i = 0; i < n; i++)
            ff[i] = func(yy, t + tau, i);
        for (int i = 0; i < n; i++)
            y[i] += tau * (a[i] + ff[i]) / 2.0;
    }
    tk = omp_get_wtime();
    deltat = tk - tn;
    printf("Массив параметров: ");
    for (int i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\nВремя вычислений: %f\n\n", deltat);
}

void RK4()
{
    const int n = 2;
    const int m = 4;
    SetConsoleOutputCP(1251);
    double t0 = 0.0, tmax = 7.0, tau = 0.01, t = t0;
    double y[n] = { 0.0, 0.5 }, yy[n] = { 0.0 }, tn, tk, deltat, R[m][n];
    tn = omp_get_wtime();
    for (double t = t0; t <= tmax; t += tau)
    {
        for (int i = 0; i < n; i++)
            R[0][i] = tau * func(y, t, i);
        for (int i = 0; i < n; i++)
            yy[i] = y[i] + 0.5 * R[0][i];
        for (int i = 0; i < n; i++)
            R[1][i] = tau * func(yy, t + 0.5 * tau, i);
        for (int i = 0; i < n; i++)
            yy[i] = y[i] + 0.5 * R[1][i];
        for (int i = 0; i < n; i++)
            R[2][i] = tau * func(yy, t + 0.5 * tau, i);
        for (int i = 0; i < n; i++)
            yy[i] = y[i] + R[2][i];
        for (int i = 0; i < n; i++)
            R[3][i] = tau * func(yy, t + tau, i);
        for (int i = 0; i < n; i++)
            y[i] += (R[0][i] + 2.0 * R[1][i] + 2.0 * R[2][i] + R[3][i]) / 6.0;
    }
    tk = omp_get_wtime();
    deltat = tk - tn;
    printf("Массив параметров: ");
    for (int i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\nВремя вычислений: %f\n\n", deltat);
}

void Euler2()
{
    const int n = 2;
    SetConsoleOutputCP(1251);
    double t0 = 0.0, tmax = 7.0, tau = 0.01, t = t0, deltah = 0.01;
    double y[n] = { 0.0, 0.5 }, tn, tk, deltat, p[n] = { 0.0 }, a[n][n] = { 0.0 }, b[n] = { 0.0 };
    double yy[n] = { y[0] + deltah, y[1]}, yy1[n] = {y[0], y[1] + deltah};
    double delta = 0.0, deltan[n] = { 0.0 };
    tn = omp_get_wtime();
    for (double t = t0; t <= tmax; t += tau)
    {
        for (int i = 0; i < n; i++)
            b[i] = -1.0 * func(y, t, i);
        for (int i = 0; i < n; i++)
            a[i][0] = (func(yy, t, i) + b[i]) / deltah;
        for (int i = 0; i < n; i++)
            a[i][1] = (func(yy1, t, i) + b[i]) / deltah;
        a[0][0] -= 1.0 / tau;
        a[1][1] -= 1.0 / tau;
        delta = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        if (!delta)
        {
            printf("delta = 0!\n");
            exit(1);
        }
        deltan[0] = b[0] * a[1][1] - b[1] * a[0][1];
        deltan[1] = b[1] * a[0][0] - b[0] * a[1][0];
        for (int i = 0; i < n; i++)
            p[i] = deltan[i] / delta;
        for (int i = 0; i < n; i++)
            y[i] += p[i];
    }
    tk = omp_get_wtime();
    deltat = tk - tn;
    printf("Массив параметров: ");
    for (int i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\nВремя вычислений: %f\n\n", deltat);
}

int main()
{
    SetConsoleOutputCP(1251);
    std::cout << "Явный метод Эйлера" << std::endl;
    Euler1();
    std::cout << "Метод Рунге-Кутта 2" << std::endl;
    RK2();
    std::cout << "Метод 'предиктор - корректор'" << std::endl;
    PredKor();
    std::cout << "Метод Рунге-Кутта 4" << std::endl;
    RK4();
    std::cout << "Неявный метод Эйлера" << std::endl;
    Euler2();
}
