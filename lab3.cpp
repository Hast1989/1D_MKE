#include <iostream>
#include <fstream>
#include <string>
#include<vector>
int numprob=1;
int numgran=1;
double a, b,ua,ub,qa,qb,Ma,Mb;
double PI = 3.14159265358979323;
int FindLeader(double** Matrix, int n, int j)
{
    int maxint = j;
    double max = 0;
    for (int i = j; i < n; i++)
    {
        if (abs(Matrix[i][j]) > max)
        {
            max = abs(Matrix[i][j]);
            maxint = i;
        }
    }
    return maxint;
}
double funcinter(double x, double y, double** coef,int** elems, double** points,int ne)
{
    double vx;
    double vy;
    double vx1;
    double vy1;
    int count1;
    double res;
    for (int i = 0; i < ne; i++)
    {
        count1 = 0;
        for (int j = 0; j < 4; j++)
        {
            vx = x - points[elems[i][j] - 1][0];
            vy = y - points[elems[i][j] - 1][1];
            vx1 = points[elems[i][(j+1)%4] - 1][0] - points[elems[i][j] - 1][0];
            vy1 = points[elems[i][(j+1)%4] - 1][1] - points[elems[i][j] - 1][1];
            res = vx1 * vy - vx * vy1;
            if (res > 0)
                count1++;
            if (count1 == 4)
                return coef[i][0] + coef[i][1] * x + coef[i][2] * y + coef[i][3] * x * y;
        }
    }
}
void SwichLines(double** Matrix, double* rightb, int i, int j)
{
    double resb = rightb[i];
    double* resline = Matrix[i];
    Matrix[i] = Matrix[j];
    Matrix[j] = resline;
    rightb[i] = rightb[j];
    rightb[j] = resb;

}
void Gauss(double** Matrix, double* rightb, int n, double* x)
{
    double d, s;
    for (int k = 0; k < n - 1; k++)
    {
        SwichLines(Matrix, rightb, FindLeader(Matrix, n, k), k);
        for (int j = k + 1; j < n; j++)
        {
            if (fabs(Matrix[j][k]) != 0)
            {
                d = Matrix[j][k] / Matrix[k][k];
                Matrix[j][k] = 0;
                for (int i = k + 1; i < n; i++)
                {
                    Matrix[j][i] = Matrix[j][i] - d * Matrix[k][i];
                }
                rightb[j] = rightb[j] - d * rightb[k];
            }
        }
    }
    for (int k = n; k > 0; k--)
    {
        d = 0;
        for (int j = k; j < n; j++)
        {
            s = Matrix[k - 1][j] * x[j];
            d = d + s;
        }
        x[k - 1] = (rightb[k - 1] - d) / Matrix[k - 1][k - 1];
    }
}
double length(double* A, double* B)
{
    return std::sqrt((A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]));
}
void vecone(double* B, double* A, double* vec)
{
    vec[0] = (A[0] - B[0]) / length(A, B);
    vec[1] = (A[1] - B[1]) / length(A, B);
    vec[2] = (A[2] - B[2]) / length(A, B);
}
void vec(double* vec,double h)
{
    vec[0] = vec[0] * h;
    vec[1] = vec[1] * h;
    vec[2] = vec[2] * h;
}
void NextPoint(double* A, double* P, double* vec)
{
    P[0] =A[0]+ vec[0];
    P[1] =A[1]+ vec[1];
    P[2] =A[2]+ vec[2];
}
void CopyPoint(double* A, double* P)
{
    P[0] = A[0];
    P[1] = A[1];
    P[2] = A[2];

}
void Polygon(int** np,double** A, double** P, int n, int m)
{
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            np[i + j * n][0] = i + 1 + j * n;
            np[i + j * n][1] = i + 1 + j * (n + 1);
            np[i + j * n][2] = i + 2 + j * (n + 1);
            np[i + j * n][3] = i + 2 + (j + 1) * (n + 1);
            np[i + j * n][4] = i + 1 + (j + 1) * (n + 1);
        }
    }
    double vec1[3];
    double vec2[3];
    double vec3[3];
    double vec4[3];
    double h1, h2,h3,h4;
    h1 = length(A[0], A[1]) / n;
    h2 = length(A[0], A[3]) / m;
    h3 = length(A[2], A[3]) / n;
    h4 = length(A[1], A[2]) / m;
    vecone(A[0], A[1], vec1);
    vecone(A[0], A[3], vec2);
    vecone(A[3], A[2], vec3);
    vecone(A[1], A[2], vec4);
    vec(vec1, h1);
    vec(vec2, h2);
    vec(vec3, h3);
    vec(vec4, h4);
    CopyPoint(A[0], P[0]);
    CopyPoint(A[1], P[n]);
    for (int j = 0; j < m; j++)
    {
        NextPoint(P[j * (n + 1)], P[(j + 1) * (n + 1)], vec2);
        NextPoint(P[n + j * (n + 1)], P[n + (j + 1) * (n + 1)], vec4);
    }
    for (int i = 0; i < n; i++)
        NextPoint(P[i + m * (n + 1)], P[i + 1 + m * (n + 1)], vec3);
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
            NextPoint(P[i+ j * (n + 1)], P[i + 1+ j * (n + 1)], vec1);
        h3 = length(P[(j+1) * (n + 1)], P[n + (j + 1) * (n + 1)]) / n;
        vecone(P[(j + 1) * (n + 1)], P[n + (j + 1) * (n + 1)], vec1);
        vec(vec1, h3);
    }
}
void Triangleprv(int** newnp3, int** np, int n, int m)
{
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            newnp3[i + j * n][0] = np[i + j * n][0];
            newnp3[i + j * n][1] = np[i + j * n][1];
            newnp3[i + j * n][2] = np[i + j * n][2];
            newnp3[i + j * n][3] = np[i + j * n][4];
            newnp3[m*n+i + j * n][0] = m*n+np[i + j * n][0];
            newnp3[m * n + i + j * n][1] = np[i + j * n][2];
            newnp3[m * n + i + j * n][2] = np[i + j * n][3];
            newnp3[m * n + i + j * n][3] = np[i + j * n][4];
        }
    }

}
void Trianglelft(int** newnp3, int** np, int n, int m)
{
    for (int j = 0; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            newnp3[i + j * n][0] = np[i + j * n][0];
            newnp3[i + j * n][1] = np[i + j * n][1];
            newnp3[i + j * n][2] = np[i + j * n][2];
            newnp3[i + j * n][3] = np[i + j * n][3];
            newnp3[m * n + i + j * n][0] = m * n + np[i + j * n][0];
            newnp3[m * n + i + j * n][1] = np[i + j * n][1];
            newnp3[m * n + i + j * n][2] = np[i + j * n][3];
            newnp3[m * n + i + j * n][3] = np[i + j * n][4];
        }
    }
}
void Triangle(int** newnp3,double** newP,int** np, double** P, int n, int m)
{
    for (int i = 0; i < n * m; i++)
    {
        newnp3[4 * i][0] = 4 * i + 1;
        newnp3[4 * i][1] = np[i][1];
        newnp3[4 * i][2] = np[i][2];
        newnp3[4 * i][3] = (n + 1) * (m + 1) +i+1;
        newnp3[4 * i+1][0] = 4 * i + 1+1;
        newnp3[4 * i+1][1] = np[i][2];
        newnp3[4 * i+1][2] = np[i][3];
        newnp3[4 * i+1][3] = (n+1) * (m+1) + i+1;
        newnp3[4 * i+2][0] = 4 * i + 1+2;
        newnp3[4 * i+2][1] = np[i][3];
        newnp3[4 * i+2][2] = np[i][4];
        newnp3[4 * i+2][3] = (n + 1) * (m + 1) + i + 1;
        newnp3[4 * i+3][0] = 4 * i + 1+3;
        newnp3[4 * i+3][1] = np[i][4];
        newnp3[4 * i+3][2] = np[i][1];
        newnp3[4 * i+3][3] = (n + 1) * (m + 1) + i+1;
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
        CopyPoint(P[i], newP[i]);
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    double u;
    for (int i = 0; i < n * m; i++)
    {
        x1 = P[np[i][1]-1][0];
        x2 = P[np[i][3]-1][0];
        x3 = P[np[i][2] - 1][0];
        x4 = P[np[i][4] - 1][0];
        y1 = P[np[i][1] - 1][1];
        y2 = P[np[i][3] - 1][1];
        y3 = P[np[i][2] - 1][1];
        y4 = P[np[i][4] - 1][1];
        z1 = P[np[i][1] - 1][2];
        z2 = P[np[i][3] - 1][2];
        z3 = P[np[i][2] - 1][2];
        z4 = P[np[i][4] - 1][2];
        u = (-x3 * y2 + x4 * y2 + x2 * y3 - x4 * y3 - x2 * y4 + x3 * y4) / (x3 * y1 - x4 * y1 - x3 * y2 + x4 * y2 - x1 * y3 + x2 * y3 + x1 * y4 - x2 * y4);
        newP[(n + 1) * (m + 1) + i][0] = u * (x1 - x2) + x2;
        newP[(n + 1) * (m + 1) + i][1] = u * (y1 - y2) + y2;
        newP[(n + 1) * (m + 1) + i][2] = u * (z1 - z2) + z2;
    }
}
void Triangle4(double** A,  int n, int m, std::ofstream& ans)
{
    double vec1[3];
    double vec2[3];
    int ne = 4*n * m;
    int np = (n + 1) * (m + 1);
    double h1, h2;
    double** P;
    double** P33;
    int** npoints;
    int** newnp33;
    P = new double* [(n + 1) * (m + 1)];
    P33 = new double* [(n + 1) * (m + 1) + m * n];
    npoints = new int* [n * m];
    newnp33 = new int* [4 * n * m];
    for (int i = 0; i < n * m; i++)
    {
        npoints[i] = new int[5];
    }
    for (int i = 0; i < 4 * n * m; i++)
    {
        newnp33[i] = new int[4];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        P[i] = new double[3];
    }
    for (int i = 0; i < (n + 1) * (m + 1) + m * n; i++)
    {
        P33[i] = new double[3];
    }
    int count;
    Polygon(npoints, A, P, n, m);
    Triangle(newnp33, P33, npoints, P, n, m);
    ans <<  ne << ' ' << np + n * m << ' ' << 1 << std::endl;
    for (int i = 0; i < 4 * n * m; i++)
    {
        ans << newnp33[i][0] <<' '<<3<<' ' << newnp33[i][1] << ' ' << newnp33[i][2] << ' ' << newnp33[i][3] << std::endl;;
    }
    for (int i = 0; i < (n + 1) * (m + 1) + n * m; i++)
    {
        ans << i + 1 << ' ' << P33[i][0] << ' ' << P33[i][1] << ' ' << P33[i][2] << std::endl;;
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        ans << i + 1 << ' ' << P[i][0] << ' ' << P[i][1] << ' ' << P[i][2] << std::endl;;
    }
    ans << 2 * n + 2 * m << std::endl;
    count = 1;
    for (int i = 0; i < n + 1; i++)
    {
        ans << count << std::endl;
        count++;
    }
    count--;
    for (int j = 1; j < m + 1; j++)
    {
        ans << count + j * (n + 1) << std::endl;
    }
    count = count + m * (n + 1);
    for (int j = 1; j < m + 1; j++)
    {
        ans << count - 1 << std::endl;
        count--;
    }
    for (int j = 1; j < m; j++)
    {
        ans << count - j * (n + 1) << std::endl;
    }
    for (int i = 0; i < n * m; i++)
    {
        delete[] npoints[i];
    }
    for (int i = 0; i < 4 * n * m; i++)
    {
        delete[] newnp33[i];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        delete[] P[i];
    }
    for (int i = 0; i < (n + 1) * (m + 1) + m * n; i++)
    {
        delete[] P33[i];
    }
    delete[] npoints;
    delete[] newnp33;
    delete[] P;
    delete[] P33;

}
void Polygon1(double** A, int n, int m, std::ofstream& ans)
{
    double vec1[3];
    double vec2[3];
    int ne = n * m;
    int np = (n + 1) * (m + 1);
    double h1, h2;
    double** P;
    int** npoints;
    P = new double* [(n + 1) * (m + 1)];
    npoints = new int* [n * m];
    for (int i = 0; i < n * m; i++)
    {
        npoints[i] = new int[5];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        P[i] = new double[3];
    }
    int count;
    Polygon(npoints, A, P, n, m);
    ans << ne << ' ' << np << ' ' << 1 << std::endl;
    for (int i = 0; i < n * m; i++)
    {
        ans << npoints[i][0] <<' '<<4<< ' ' << npoints[i][1] << ' ' << npoints[i][2] << ' ' << npoints[i][3] << ' ' << npoints[i][4] << std::endl;;
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        ans << i + 1 << ' ' << P[i][0] << ' ' << P[i][1] << ' ' << P[i][2] << std::endl;
    }
    ans << 2 * n + 2 * m << std::endl;
    count = 1;
    for (int i = 0; i < n + 1; i++)
    {
        ans << count << std::endl;
        count++;
    }
    count--;
    for (int j = 1; j < m + 1; j++)
    {
        ans << count + j * (n + 1) << std::endl;
    }
    count = count + m * (n + 1);
    for (int j = 1; j < m + 1; j++)
    {
        ans << count - 1 << std::endl;
        count--;
    }
    for (int j = 1; j < m; j++)
    {
        ans << count - j * (n + 1) << std::endl;
    }
    for (int i = 0; i < n * m; i++)
    {
        delete[] npoints[i];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        delete[] P[i];
    }
    delete[] npoints;
    delete[] P;
}
void Triangleprv1(double** A, int n, int m, std::ofstream& ans)
{
    double vec1[3];
    double vec2[3];
    int ne =2* n * m;
    int np = (n + 1) * (m + 1);
    double h1, h2;
    double** P;
    int** npoints;
    int** newnp3;
    P = new double* [(n + 1) * (m + 1)];
    npoints = new int* [n * m];
    newnp3 = new int* [2 * n * m];
    for (int i = 0; i < n * m; i++)
    {
        npoints[i] = new int[5];
    }
    for (int i = 0; i < 2 * n * m; i++)
    {
        newnp3[i] = new int[4];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        P[i] = new double[3];
    }
    Polygon(npoints, A, P, n, m);
    Triangleprv(newnp3, npoints, n, m);
    ans << ne << ' ' << np << ' ' << 1 << std::endl;
    for (int i = 0; i <2* n * m; i++)
    {
        ans << newnp3[i][0] <<' ' << 3 << ' ' << newnp3[i][1] << ' ' << newnp3[i][2] << ' ' << newnp3[i][3] << std::endl;;
    }
    int count;
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        ans << i + 1 << ' ' << P[i][0] << ' ' << P[i][1] << ' ' << P[i][2] << std::endl;;
    }
    ans << 2 * n + 2 * m << std::endl;
    count = 1;
    for (int i = 0; i < n + 1; i++)
    {
        ans << count << std::endl;
        count++;
    }
    count--;
    for (int j = 1; j < m + 1; j++)
    {
        ans << count + j * (n + 1) << std::endl;
    }
    count = count + m * (n + 1);
    for (int j = 1; j < m + 1; j++)
    {
        ans << count - 1 << std::endl;
        count--;
    }
    for (int j = 1; j < m; j++)
    {
        ans << count - j * (n + 1) << std::endl;
    }
    for (int i = 0; i < n * m; i++)
    {
        delete[] npoints[i];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        delete[] P[i];
    }
    delete[] npoints;
    delete[] P;
}
void Trianglelft1(double** A, int n, int m, std::ofstream& ans)
{
    double vec1[3];
    double vec2[3];
    int ne =2* n * m;
    int np = (n + 1) * (m + 1);
    double h1, h2;
    double** P;
    int** npoints;
    int** newnp3;
    P = new double* [(n + 1) * (m + 1)];
    npoints = new int* [n * m];
    newnp3 = new int* [2 * n * m];
    for (int i = 0; i < n * m; i++)
    {
        npoints[i] = new int[5];
    }
    for (int i = 0; i < 2 * n * m; i++)
    {
        newnp3[i] = new int[4];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        P[i] = new double[3];
    }
    Polygon(npoints, A, P, n, m);
    Trianglelft(newnp3, npoints, n, m);
    ans << ne << ' ' << np << ' ' << 1 << std::endl;
    for (int i = 0; i < 2 * n * m; i++)
    {
        ans << newnp3[i][0] << ' ' << 3<< ' ' << newnp3[i][1] << ' ' << newnp3[i][2] << ' ' << newnp3[i][3] << std::endl;;
    }
    int count;
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        ans << i + 1 << ' ' << P[i][0] << ' ' << P[i][1] << ' ' << P[i][2] << std::endl;;
    }
    ans << 2 * n + 2 * m << std::endl;
    count = 1;
    for (int i = 0; i < n + 1; i++)
    {
        ans << count << std::endl;
        count++;
    }
    count--;
    for (int j = 1; j < m + 1; j++)
    {
        ans << count + j * (n + 1) << std::endl;
    }
    count = count + m * (n + 1);
    for (int j = 1; j < m + 1; j++)
    {
        ans << count - 1 << std::endl;
        count--;
    }
    for (int j = 1; j < m; j++)
    {
        ans << count - j * (n + 1) << std::endl;
    }
    for (int i = 0; i < n * m; i++)
    {
        delete[] npoints[i];
    }
    for (int i = 0; i < (n + 1) * (m + 1); i++)
    {
        delete[] P[i];
    }
    delete[] npoints;
    delete[] P;
}
void Linear(double** A, int n, std::ofstream& ans)
{
    double h;
    h = length(A[0], A[1])/n;
    int ne = n;
    int np = n + 1;
    double** P;
    int** npoints;
    double vec1[3];
    P = new double* [(n + 1)];
    npoints = new int* [n];
    for (int i = 0; i < n; i++)
    {
        npoints[i] = new int[4];
    }
    for (int i = 0; i < (n + 1); i++)
    {
        P[i] = new double[3];
    }
    vecone(A[0], A[1], vec1);
    vec(vec1, h);
    CopyPoint(A[0], P[0]);
    for (int i = 0; i < n; i++)
    {
        NextPoint(P[i], P[i + 1], vec1);
        npoints[i][0] = i + 1;
        npoints[i][1] = 2;
        npoints[i][2] = i + 1;
        npoints[i][3] = i + 2;
    }
    ans << ne << ' ' << np << ' ' << 1 << std::endl;
    for (int i = 0; i < n; i++)
    {
        ans << npoints[i][0] << ' ' << npoints[i][1] << ' ' << npoints[i][2] << ' ' << npoints[i][3] << std::endl;;
    }
    for (int i = 0; i < (n + 1); i++)
    {
        ans << i + 1 << ' ' << P[i][0] << ' ' << P[i][1] << ' ' << P[i][2] << std::endl;;
    }
    ans << 1 <<' ' <<1<< std::endl;
    ans << 1 << std::endl;
    ans << 5 << std::endl;
    for (int i = 0; i < n; i++)
    {
        delete[] npoints[i];
    }
    for (int i = 0; i < (n + 1); i++)
    {
        delete[] P[i];
    }
    delete[] npoints;
    delete[] P;

}
void SqLinear(double** A, int n, std::ofstream& ans)
{
    double h;
    h = length(A[0], A[1]) / n;
    h = h / 2;
    int ne = n;
    int np = 2*n + 1;
    double** P;
    int** npoints;
    double vec1[3];
    P = new double* [2*n + 1];
    npoints = new int* [n];
    for (int i = 0; i < n; i++)
    {
        npoints[i] = new int[5];
    }
    for (int i = 0; i < 2*n+1; i++)
    {
        P[i] = new double[3];
    }
    vecone(A[0], A[1], vec1);
    vec(vec1, h);
    CopyPoint(A[0], P[0]);
    for (int i = 0; i < 2*n; i++)
    {
        NextPoint(P[i], P[i + 1], vec1);
        if (i % 2 == 0)
        {
            npoints[i/2][0] = i/2 + 1;
            npoints[i/2][1] = 2;
            npoints[i/2][2] = i + 1;
            npoints[i/2][3] = i + 2;
            npoints[i/2][4] = i + 3;
        }
    }
    ans << ne << ' ' << np << ' ' << 1 << std::endl;
    for (int i = 0; i < n; i++)
    {
        ans << npoints[i][0] << ' ' <<  npoints[i][1] << ' ' << npoints[i][2] << ' ' << npoints[i][3] << ' ' << npoints[i][4] << std::endl;;
    }
    for (int i = 0; i < 2*n + 1; i++)
    {
        ans << i + 1 << ' ' << P[i][0] << ' ' << P[i][1] << ' ' << P[i][2] << std::endl;
    }
    ans << 1 << ' ' << 1 << std::endl;
    ans << 1 << std::endl;
    ans << 5 << std::endl;
    for (int i = 0; i < n; i++)
    {
        delete[] npoints[i];
    }
    for (int i = 0; i < (n + 1); i++)
    {
        delete[] P[i];
    }
    delete[] npoints;
    delete[] P;

}
double func(double x, double y, double z)
{
    return 10 + 5 * y * std::sin(2 * PI * x);
}
void Interpabcdel(int* elem, double** points, double* f,double* coef)
{
    double** A;
    A = new double* [4];
    double* PR;
    PR = new double[4];
    for (int i = 0; i < 4; i++)
    {
        A[i] = new double[4];
        A[i][0] = 1;
        A[i][1] = points[elem[i] - 1][0];
        A[i][2] = points[elem[i] - 1][1];
        A[i][3] = points[elem[i] - 1][0] * points[elem[i] - 1][1];
        PR[i] = f[elem[i] - 1];
    }
    Gauss(A, PR, 4, coef);
    delete[] PR;
    for (int i = 0; i < 4; i++)
    {
        delete[] A[i];
    }
    delete[] A;
}
double Intrel(int* elem, double** points, double* f)
{
    double v1x = points[elem[1] - 1][0] - points[elem[0] - 1][0];
    double v1y = points[elem[1] - 1][1] - points[elem[0] - 1][1];
    double v2x = points[elem[3] - 1][0] - points[elem[0] - 1][0];
    double v2y = points[elem[3] - 1][1] - points[elem[0] - 1][1];
    double v3x = points[elem[1] - 1][0] - points[elem[2] - 1][0];
    double v3y = points[elem[1] - 1][1] - points[elem[2] - 1][1];
    double v4x = points[elem[3] - 1][0] - points[elem[2] - 1][0];
    double v4y = points[elem[3] - 1][1] - points[elem[2] - 1][1];
    double S = std::fabs(v1x * v2y - v1y * v2x) +std::fabs(v3x * v4y - v3y * v4x);
    double zn = (f[elem[0]-1] + f[elem[1]-1] + f[elem[2]-1] + f[elem[3]-1]) / 4;
    return S * zn*0.5;
}
std::vector<int> findneig(int ne,int** elems, int nump)
{
    std::vector<int> res;
    for (int i = 0; i < ne; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            if (nump == elems[i][j])
            {
                res.push_back(i + 1);
            }
        }
    }
    return res;
}
void test(std::string in,std::string out)
{
    int res;
    int meth;
    std::ifstream data(in);
    std::ofstream file;
    file.open(out);
    data >> res;
    if (res == 4)
    {
        int n;
        int m;
        double** A;
        A = new double* [4];
        for (int i = 0; i < 4; i++)
        {
            A[i] = new double[3];
        }
        data >> A[0][0];
        data >> A[0][1];
        data >> A[0][2];
        data >> A[1][0];
        data >> A[1][1];
        data >> A[1][2];
        data >> A[2][0];
        data >> A[2][1];
        data >> A[2][2];
        data >> A[3][0];
        data >> A[3][1];
        data >> A[3][2];
        data >> n;
        data >> m;
        data >> meth;
        if (meth == 1)
            Polygon1(A, n, m, file);
        if (meth == 2)
            Trianglelft1(A, n, m, file);
        if (meth == 3)
            Triangleprv1(A, n, m, file);
        if (meth == 4)
            Triangle4(A, n, m, file);
        file.close();
    }
    if (res == 2)
    {
        int n;
        double** A;
        A = new double* [2];
        for (int i = 0; i < 2; i++)
        {
            A[i] = new double[3];
        }
        data >> A[0][0];
        data >> A[0][1];
        data >> A[0][2];
        data >> A[1][0];
        data >> A[1][1];
        data >> A[1][2];
        data >> n;
        data >> meth;
        if (meth == 1)
            Linear(A, n, file);
        if (meth == 2)
            SqLinear(A, n, file);
    }
    data.close();
    file.close();
}
void testfunc(std::string in, std::string out1, std::string out2, std::string out3, std::string out4, std::string out5)
{
    int np, ne;
    int res;
    std::ifstream data(in);
    std::ofstream file1;
    std::ofstream file2;
    std::ofstream file3;
    std::ofstream file4;
    std::ofstream file5;
    data >> ne;
    data >> np;
    data >> res;
    int** elems;
    double** coef;
    double** points;
    double* f;
    double* ife;
    double* ifp;
    double** pfe;
    double** pfp;
    elems = new int* [ne];
    coef = new double* [ne];
    points = new double* [np];
    pfe = new double* [ne];
    pfp = new double* [np];
    f= new double[np];
    ife = new double[ne];
    ifp = new double[np];
    in = "out.txt";
    out1 = "fzn.txt";
    out2 = "prel.txt";
    out3 = "iel.txt";
    out4 = "pruz.txt";
    out5 = "iuz.txt";
    for (int i = 0; i < ne; i++)
    {
        data >> res;
        data >> res;
        pfe[i] = new double[2];
        coef[i] = new double[res];
        elems[i] = new int[res];
        for (int j = 0; j < res; j++)
        {
            data >> elems[i][j];
        }
    }
    for (int i = 0; i < np; i++)
    {
        data >> res;
        pfp[i] = new double[2];
        points[i] = new double[3];
        for (int j = 0; j < 3; j++)
        {
            data >> points[i][j];
        }
    }
    data.close();
    file1.open(out1);
    file2.open(out2);
    file3.open(out3);
    file4.open(out4);
    file5.open(out5);
    file1 << np << std::endl;
    for (int i = 0; i < np; i++)
    {
        f[i] = func(points[i][0], points[i][1], points[i][2]);
        file1 << points[i][0] << ' ' << points[i][1] << ' ' << points[i][2] << ' ' << f[i] << std::endl;
    }
    
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    double u;
    file1 << ne << std::endl;
    for (int i = 0; i < ne; i++)
    {
        Interpabcdel(elems[i], points, f, coef[i]);
        x1 = points[elems[i][0]-1][0];
        x2 = points[elems[i][1]-1][0];
        x3 = points[elems[i][2]-1][0];
        x4 = points[elems[i][3]-1][0];
        y1 = points[elems[i][0]-1][1];
        y2 = points[elems[i][1]-1][1];
        y3 = points[elems[i][2]-1][1];
        y4 = points[elems[i][3]-1][1];
        u = (-x3 * y2 + x4 * y2 + x2 * y3 - x4 * y3 - x2 * y4 + x3 * y4) / (x3 * y1 - x4 * y1 - x3 * y2 + x4 * y2 - x1 * y3 + x2 * y3 + x1 * y4 - x2 * y4);
        pfe[i][0] = coef[i][1] + coef[i][3] * (u * (y1 - y2) + y2);
        pfe[i][1] = coef[i][2] + coef[i][3] * (u * (x1 - x2) + x2);
        file1 << coef[i][0] << ' ' << coef[i][1] << ' ' << coef[i][2] << ' ' << coef[i][3] << std::endl;
        file2 << "( " << coef[i][1] << " + " << coef[i][3] << " * y ) dx" << "+ ( " << coef[i][2] << " + " << coef[i][3] << " * x ) dy" << std::endl;
        file2 << pfe[i][0]<<' '<< pfe[i][1] << std::endl;
        ife[i] = Intrel(elems[i], points, f);
        file3 << ife[i] << std::endl;
    }
    std::cout << funcinter(5, 5, coef, elems, points, ne) << std::endl;
    for (int i = 0; i < np; i++)
    {
        std::vector<int> res;
        res = findneig(ne, elems, i + 1);
        double sum1 = 0;
        double sum2 = 0;
        double sumI = 0;
        for (int j = 0; j < res.size(); j++)
        {
            sum1 += pfe[res[j] - 1][0];
            sum2 += pfe[res[j] - 1][1];
            sumI += ife[res[j]-1];
        }
        sum1 = sum1 / res.size();
        sum2 = sum2 / res.size();
        sumI = sumI / res.size();
        file4 << sum1 << ' ' << sum2 << std::endl;
        file5 << sumI << std::endl;
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();
    for (int i = 0; i < ne; i++)
    {
        delete[] pfe[i];
        delete[]  coef[i];
        delete[] elems[i];
    }
    points = new double* [np];
    for (int i = 0; i < np; i++)
    {
        delete[] pfp[i];
        delete[] points[i];
    }
    delete[] elems;
    delete[] points;
    delete[] coef;
    delete[] ife;
    delete[] ifp;
    delete[] f;
    delete[] pfp;
    delete[] ifp;
}
double funclammda(double x)
{
    if (numprob == 1)
    {
        return 10+20 * x * (x - 1);
    }
    if (numprob == 2)
    {
        return 10;
    }
}
double funcQ(double x)
{
    if (numprob == 1)
    {
        return 10000*(x-0.5);
    }
    if (numprob == 2)
    {
        return -10000 * (x - 0.5);
    }
}
double funcalphaa()
{
    if (numprob == 1)
    {
        return 100000;
    }
    if (numprob == 2)
    {
        return 100;
    }
}
double funcalphab()
{
    if (numprob == 1)
    {
        return 100;
    }
    if (numprob == 2)
    {
        return 100;
    }
}
double uinfa()
{
    if (numprob == 1)
    {
        return 200;
    }
    if (numprob == 2)
    {
        return 300;
    }
}
double uinfb()
{
    if (numprob == 1)
    {
        return 300;
    }
    if (numprob == 2)
    {
        return 300;
    }
}
void testel(std::string in, std::string out)
{
    numprob = 2;
    numgran = 1;
    std::ifstream data(in);
    std::ofstream file(out);
    int ne;
    int np;
    double* k;
    double* uz;
    double* u;
    double** KG;
    double* B;
    data >> ne;
    data >> np;
    uz = new double[np];
    k = new double[4 * np];
    u= new double[np];
    B = new double[np];
    KG = new double* [ne];
    data >> k[0];
    for (int i = 0; i < 4 * ne; i++)
        data >> k[i];
    for (int i = 0; i < 4 * np; i++)
        data >> k[i];
    for (int i = 0; i < np; i++)
    {
        uz[i] = k[4 * i + 1];
        KG[i] = new double[np];
    }
    for (int i = 0; i < np; i++)
        for (int j = 0; j < np; j++)
            KG[i][j] = 0.;
    KG[0][0] = ((funclammda(uz[1]) + funclammda(uz[0])) / 2) / (uz[1] - uz[0]);
    KG[0][1] = -((funclammda(uz[1]) + funclammda(uz[0])) / 2) / (uz[1] - uz[0]);
    B[0] =funcQ((uz[1] + uz[0]) / 2)* (uz[1] - uz[0])/2;
    for (int i = 1; i < np-1; i++)
    {
        KG[i][i-1] = -((funclammda(uz[i])+ funclammda(uz[i-1]))/2) / (uz[i] - uz[i-1]);
        KG[i][i] = ((funclammda(uz[i]) + funclammda(uz[i - 1])) / 2) / (uz[i] - uz[i-1]) + ((funclammda(uz[i+1]) + funclammda(uz[i])) / 2) / (uz[i + 1] - uz[i]);
        KG[i][i+1] = -((funclammda(uz[i + 1]) + funclammda(uz[i])) / 2) / (uz[i + 1] - uz[i]);
        B[i] = funcQ((uz[i + 1] + uz[i]) / 2)* (uz[i + 1] - uz[i])/2+  funcQ((uz[i] + uz[i-1]) / 2) * (uz[i] - uz[i-1])/2;
    }
    KG[np-1][np-2] = -((funclammda(uz[np-1]) + funclammda(uz[np-2])) / 2) / (uz[np-1] - uz[np-2]);
    KG[np - 1][np - 1] = ((funclammda(uz[np - 1]) + funclammda(uz[np - 2])) / 2) / (uz[np - 1] - uz[np - 2]);
    B[np-1] =  funcQ((uz[np-1] + uz[np-2]) / 2) * (uz[np-1] - uz[np-2])/2;
    if (numgran == 1)
    {
        KG[0][0] -= Ma * (uz[1] - uz[0]) / 2;
        KG[np-1][np-1] -= Mb * (uz[np-1] - uz[np-2]) / 2;
        B[0] -= Ma * ua * (uz[1] - uz[0]) / 2;
        B[np - 1]-= Mb * ub * (uz[np - 1] - uz[np - 2]) / 2;
    }
    if (numgran == 2)
    {
        B[0] += qa;
        B[np - 1] -= qb;
    }
    if (numgran == 3)
    {
        B[0] += funcalphaa()*uinfa() ;
        B[np - 1] += funcalphab()*uinfb() ;
        KG[0][0]+= funcalphaa() ;
        KG[np - 1][np - 1]+= funcalphab() ;
    }
    Gauss(KG, B, np, u);
    for (int i = 0; i < np; i++)
        file << uz[i] << ' ' << u[i] << std::endl;
    
}       
int main()
{
    Ma = 9000000;
    Mb = 9000000;
    ua = 7;
    ub = 10;
    qa = 5;
    qb = 4;
    std::string in = "in.txt";
    std::string out = "out.txt";
    test(in,out);
    //std::string outf = "outf.txt";
    //testel(out, outf);
    std::cout << "Hello World!\n";
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
