#include "mainwindow.h"
#include "./ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{            
    ui->setupUi(this);  
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::plotchart(){
    m_axisX->setTitleText("x坐标");
    m_axisX->setRange(X_MIN , X_MAX);
    m_axisY->setTitleText("y坐标");
    m_axisY->setRange(Y_MIN , Y_MAX);

    m_series->setPointsVisible(true);
    m_series->setName("迭代值");

    m_chart->addSeries(m_series);
    m_chart->setTitle("迭代解的变化");
    m_chart->addAxis(m_axisX,Qt::AlignBottom);
    m_chart->addAxis(m_axisY,Qt::AlignLeft);
    m_chart->legend()->setAlignment(Qt::AlignBottom);
    m_chart->setAnimationOptions(QChart::SeriesAnimations);

    ui->graphicsView->setChart(m_chart);
    ui->graphicsView->setRenderHint(QPainter::Antialiasing); //抗锯齿
}

void MainWindow::updateMessage(QString& message) {
    ui->textBrowser->append(message);
}

void MainWindow::updateMessage(QString&& message) {
    ui->textBrowser->append(message);
}

void MainWindow::updateChart(const double& x, const double& y){
    //调整X和Y轴范围
    if(m_series->count() == 0){
        X_MIN = x; X_MAX = x;
        Y_MIN = y; Y_MAX = y;
    } else {
        if(x < X_MIN)
            X_MIN = x;
        if(x > X_MAX)
            X_MAX = x;
        if(y < Y_MIN)
            Y_MIN = y;
        if(y > Y_MAX)
            Y_MAX = y;
    }
    m_series->append(x, y);
}




