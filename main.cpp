#include "mainwindow.h"
#include <QApplication>
#include "LaLib.h"

using std::vector;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

/*********************一些样板矩阵********************/
//    //二阶 谱半径小于1   对称
//    vector<vector<double>> A = {{0.5, 0.3},
//                                {0.3, 0.5}};
//    vector<double> b = {8, 8};
//    vector<double> x = {0, 0};

    //二阶 谱半径大于1    对称
    vector<vector<double>> A = {{3, 2},
                                {2, 6}};
    vector<double> b = {2, -8};
    vector<double> x = {0, 0};


//    //二阶 谱半径小于1   不对称
//    vector<vector<double>> A = {{0.5, 0.3},
//                                {0.2, 0.5}};
//    vector<double> b = {8, 7};
//    vector<double> x = {0, 0};


//    //二阶 谱半径大于1    不对称
//    vector<vector<double>> A = {{3, 5},
//                                {1, 6}};
//    vector<double> b = {2, -8};
//    vector<double> x = {0, 0};


//    //三阶 不对称
//    vector<vector<double>> A = {{3, 5, 3},
//                                {1, 6, 5},
//                                {2, 4, 6}};
//    vector<double> b = {2, -8, 1};
//    vector<double> x = {0, 0, 0};


//    //五阶 不对称
//    vector<vector<double>> A = {{2.175,  -0.5917, 0,       0,       0     },
//                                {-0.7 ,  1.075,   -0.425,  0,       0     },
//                                {0.025,  -0.675,  1.075,   -0.425,  0     },
//                                {0,      0.025,   -0.675,  1.075,   -0.425},
//                                {0,      0,       0.025,   -0.8167, 1.925 }};
//    vector<double> b = {1.583, -0.05, 0, 0, 0};
//    vector<double> x = {0, 0, 0, 0, 0};

/*********************线性代数库的使用********************/
//    LinearAlgebra::jacobian(A, x, b, w);    //当谱半径大于1时雅可比迭代很可能不收敛（二阶矩阵仍有可能收敛）
//    LinearAlgebra::GaussSeidel(A, x, b, w);
//    LinearAlgebra::jacobianSOR(A, x, b, 1.5, w);   //超松弛
//    LinearAlgebra::jacobianSOR(A, x, b, 0.7, w);   //欠松弛
//    LinearAlgebra::GaussSeidelSOR(A, x, b, 1.5, w); //超松弛
//    LinearAlgebra::GaussSeidelSOR(A, x, b, 0.7, w); //欠松弛
    LinearAlgebra::GradientDescent(A, x, b, w);     //当矩阵阶数小时非对称矩阵，梯度下降法仍有可能得到解
//    LinearAlgebra::CG(A, x, b, w);  //当矩阵非对称时CG体现在最后阶段很难收敛下去
//    LinearAlgebra::BICGStab(A, x, b, w);

/**********调用画图函数（只有矩阵为2阶时才有数据）**************/
    w.plotchart();
    return a.exec();
}
