#ifndef LALIB_H
#define LALIB_H
#include <iostream>
#include <cmath>
#include <QtCharts>
#include <vector>
#include "mainwindow.h"

using std::vector;

/********** 下面是生成QString的辅助函数 **********/
QString vectorToQString(const vector<double>& v){
    QString result = "[     ";
    for(int i = 0; i < v.size(); i++){
        //https://stackoverflow.com/questions/65199469/format-qstring-with-fixed-number-of-digits
        result += QString("%1").arg(v[i], 5, 'f', 7) + "     ";
    }
    result += "]";
    return result;
}

namespace LinearAlgebra{
//iteration parameter
const int ITERATION_LIMIT = 100;
const double rtol = 1e-5;

/********** 下面是一些线性代数的辅助函数 **********/
QString EigValForTwoOrderMatrix(const vector<vector<double>>& A){
    QString result;
    double m = (A[0][0] + A[1][1]) / 2;
    double p = A[0][0] * A[1][1] - A[0][1] * A[1][0];   
    double m2_minus_p = m * m - p;
    if(m2_minus_p < 0){
        result = "No Eigen Value";
    } else if(m2_minus_p == 0) {
        result = "one Eigen Value: " + QString::number(m);
    } else {
        result = "two Eigen Value      Eigen Value1:  " + QString::number(m - sqrt(m2_minus_p)) +
                 "       Eigen Value2:  " + QString::number(m + sqrt(m2_minus_p));
    }
    return result;
}

vector<double> linearCombine(const double a1, const vector<double>& v1, const double a2, const vector<double>& v2) {
    //输入检测机制待插入
    vector<double> result;
    for(int i = 0; i < v1.size(); i++){
        result.push_back(a1 * v1[i] + a2 * v2[i]);
    }
    return result;
}

double dotProducts(const vector<double>& v1, const vector<double>& v2) {
    //输入检测机制待插入
    double result = 0;
    for(int i = 0; i < v1.size(); i++){
        result += v1[i] * v2[i];
    }
    return result;
}

vector<double> matrixMultiplication(const vector<vector<double>>& A, const vector<double>& x) {
    //输入检测机制待插入
    vector<double> result;
    for(int i = 0; i < x.size(); i++){
        double temp = 0;
        for(int j = 0; j < x.size(); j++){
            temp += A[i][j] * x[j];
        }
        result.push_back(temp);
    }
    return result;
}

bool ifconverge(const vector<double>& Newx, const vector<double>& x) {
    for (int i = 0; i < x.size(); i++){
        if (fabs(Newx[i] - x[i]) >= rtol* fabs(Newx[i])){
            return false;
        }
    }
    return true;
}

/********** 下面是一些线性代数的迭代方法 **********/
/**** 以后需要改为bool类型并加判断函数判断是否发散***/
void jacobian(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值,并更新图形信息
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }
    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    //1.初始化迭代步数和两个求和数
    //拷贝构造，创建新的vector<double> Newx
    int step = 0;
    double sum1 = 0, sum2 = 0;
    vector<double> Newx(x);
    while (true){
        //2.更新迭代步数并判断有没有达到最大迭代步数
        step++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //3.计算Newx
        for (int i = 0; i < x.size(); i++){
            for (int j = 0; j < i ; j++){
                sum1 += A[i][j] * x[j];
            }
            for (int j = i+1; j < x.size(); j++){
                sum2 += A[i][j] * x[j];
            }
            Newx[i] = (1.0 / A[i][i]) * (b[i] - sum1 - sum2);
            //4.重制sum
            sum1 = 0; sum2 = 0;
        }
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(Newx);
        w.updateMessage(message);
        //图形信息，如果阶数为2，则将数据加入数据序列
        if(x.size() == 2){
            w.updateChart(x[0], x[1]);
        }
        //5.判断是否收敛
        if(ifconverge(Newx, x))
            break;
        //6.更新迭代值
        x = Newx;
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}

void GaussSeidel(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值,并更新图形信息
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }
    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    //1.初始化迭代步数和两个求和数
    //拷贝构造，创建新的vector<double> Newx
    int step = 0;
    double sum1 = 0, sum2 = 0;
    vector<double> Newx(x);
    while (true){
        //2.更新迭代步数并判断有没有达到最大迭代步数
        step++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //3.计算Newx
        for (int i = 0; i < x.size(); i++){
            for (int j = 0; j < i ; j++){
                sum1 += A[i][j] * Newx[j];
            }
            for (int j = i+1; j < x.size(); j++){
                sum2 += A[i][j] * x[j];
            }
            Newx[i] = (1.0 / A[i][i]) * (b[i] - sum1 - sum2);
            //4.重制sum
            sum1 = 0; sum2 = 0;
        }
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(Newx);
        w.updateMessage(message);
        //图形信息，如果阶数为2，则将数据加入数据序列
        if(x.size() == 2){
            w.updateChart(x[0], x[1]);
        }
        //5.判断是否收敛
        if(ifconverge(Newx, x))
            break;
        //6.更新迭代值
        x = Newx;
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}

void jacobianSOR(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, const double& relaxtion, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值,并更新图形信息
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }
    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    //1.初始化迭代步数和两个求和数
    //拷贝构造，创建新的vector<double> Newx
    int step = 0;
    double sum1 = 0, sum2 = 0;
    vector<double> Newx(x);
    while (true){
        //2.更新迭代步数并判断有没有达到最大迭代步数
        step++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //3.计算Newx
        for (int i = 0; i < x.size(); i++){
            for (int j = 0; j < i ; j++){
                sum1 += A[i][j] * x[j];
            }
            for (int j = i+1; j < x.size(); j++){
                sum2 += A[i][j] * x[j];
            }
            Newx[i] = (1.0 - relaxtion) * x[i] + (relaxtion / A[i][i]) * (b[i] - sum1 - sum2);
            //4.重制sum
            sum1 = 0; sum2 = 0;
        }
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(Newx);
        w.updateMessage(message);
        //图形信息，如果阶数为2，则将数据加入数据序列
        if(x.size() == 2){
            w.updateChart(x[0], x[1]);
        }
        //5.判断是否收敛
        if(ifconverge(Newx, x))
            break;
        //6.更新迭代值
        x = Newx;
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}

void GaussSeidelSOR(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, const double& relaxtion, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值,并更新图形信息
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }
    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    //1.初始化迭代步数和两个求和数
    //拷贝构造，创建新的vector<double> Newx
    int step = 0;
    double sum1 = 0, sum2 = 0;
    vector<double> Newx(x);
    while (true){
        //2.更新迭代步数并判断有没有达到最大迭代步数
        step++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //3.计算Newx
        for (int i = 0; i < x.size(); i++){
            for (int j = 0; j < i ; j++){
                sum1 += A[i][j] * Newx[j];
            }
            for (int j = i+1; j < x.size(); j++){
                sum2 += A[i][j] * x[j];
            }
            Newx[i] = (1.0 - relaxtion) * x[i] + (relaxtion / A[i][i]) * (b[i] - sum1 - sum2);
            //4.重制sum
            sum1 = 0; sum2 = 0;
        }
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(Newx);
        w.updateMessage(message);
        //图形信息，如果阶数为2，则将数据加入数据序列
        if(x.size() == 2){
            w.updateChart(x[0], x[1]);
        }
        //5.判断是否收敛
        if(ifconverge(Newx, x))
            break;
        //6.更新迭代值
        x = Newx;
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}

void GradientDescent(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值,并更新图形信息
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }

    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    int step = 0;
    //1.
    vector<double> r = linearCombine(1, b, -1, matrixMultiplication(A, x));
    while(true){
        //2.
        step ++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //3.
        vector<double> Ar = matrixMultiplication(A, r);
        //4.
        double rTr = dotProducts(r, r);
        //除数不能为0
        if(rTr == 0)
            break;
        //5.
        double Gamma = rTr / dotProducts(r, Ar);
        //6.
        x = linearCombine(1, x, Gamma, r);
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(x);
        w.updateMessage(message);
        //如果阶数为2，则将数据加入数据序列
        if(x.size() == 2)
            w.updateChart(x[0], x[1]);
        //7.如果rTr够小
        if(sqrt(rTr) < rtol){
            break;
        }
        //8.更新残差
        r = linearCombine(1, r, -Gamma, Ar);
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}

void CG(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }
    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    int step = 0;
    //1.
    vector<double> r = linearCombine(1, b, -1, matrixMultiplication(A, x));
    //2.
    vector<double> p(r);
    //3.
    double rTrOld = dotProducts(r, r);
    while(true){
        //4.
        step ++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //5.
        vector<double> Ap = matrixMultiplication(A, p);
        //6.
        double alpha = rTrOld / dotProducts(p, Ap);
        //7.
        x = linearCombine(1, x, alpha, p);
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(x);
        w.updateMessage(message);
        //如果阶数为2，则将数据加入数据序列
        if(x.size() == 2)
            w.updateChart(x[0], x[1]);

        //8.清除舍入误差
        if(x.size() <= 10000 && step % 50 == 0){
            r = linearCombine(1, b, -1, matrixMultiplication(A, x));
        } else if(x.size() > 10000 && step % static_cast<int>(sqrt(x.size())) == 0){
            r = linearCombine(1, b, -1, matrixMultiplication(A, x));
        } else {
            r = linearCombine(1, r, -alpha, Ap);
        }
        //9.
        double rTrNew = dotProducts(r, r);
        if(sqrt(rTrNew) < rtol)
            break;
        //10.
        p = linearCombine(1, r, (rTrNew / rTrOld), p);
        //11.
        rTrOld = rTrNew;
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}

void BICGStab(const vector<vector<double>>& A, vector<double>& x, const vector<double>& b, MainWindow& w){
    QString message;
    //如果矩阵阶数为2，则显示特征值
    if(x.size() == 2){
        message = EigValForTwoOrderMatrix(A);
        w.updateMessage(message);
        w.updateChart(x[0], x[1]);
    }
    //显示初始迭代值
    message = "step0   " + vectorToQString(x);
    w.updateMessage(message);

    int step = 0;
    //1.
    vector<double> r = linearCombine(1, b, -1, matrixMultiplication(A, x));
    //2.
    vector<double> r0(r);
    //3.
    double bSqnorm = dotProducts(b, b);
    //如果b向量的平方范数为0，说明方程右侧全为0，方程无解
    if(bSqnorm == 0){
        x.assign(x.size(), 0);
        return;
    }
    //4.
    double rho = 1, alpha =1, omega =1;
    //5.
    vector<double> v(x.size(), 0);
    vector<double> p(x.size(), 0);
    //6.
    while(dotProducts(r, r) > rtol * rtol* bSqnorm){
        step ++;
        if (step >= ITERATION_LIMIT){
            w.updateMessage("Not converge!!!!");
            return;
        }
        //7.
        double rhoOld = rho;
        //8.
        rho = dotProducts(r0, r);
        //9.
        double beta = (rho / rhoOld) * (alpha / omega);
        //10.
        p = linearCombine(1, r0, beta, linearCombine(1, p, -omega, v));
        //11.
        v = matrixMultiplication(A, p);
        //12.
        alpha = rho / dotProducts(r0, v);
        //13.
        vector<double> h = linearCombine(1, x, alpha, p);
        //14.
        vector<double> s = linearCombine(1, r, -alpha, v);
        //15.
        vector<double> t = matrixMultiplication(A, s);
        //16.防止除数 = 0
        if(dotProducts(t, t) == 0){
            omega = 0;
        } else {
            omega = dotProducts(t, s) / dotProducts(t, t);
        }
        //17.
        x = linearCombine(1, h, omega, s);
        //更新信息
        message = "step" +  QString::number(step) + "   " + vectorToQString(x);
        w.updateMessage(message);
        //如果阶数为2，则将数据加入数据序列
        if(x.size() == 2)
            w.updateChart(x[0], x[1]);
        //18.
        r = linearCombine(1, s, -omega, t);
    }
    //显示最终迭代误差
    vector<double> error = linearCombine(1, matrixMultiplication(A, x), -1, b);
    message = "Error:     " +  vectorToQString(error);
    w.updateMessage(message);
}
}
#endif // LALIB_H
