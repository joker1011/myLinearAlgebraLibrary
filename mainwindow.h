#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <vector>
#include <QtCharts>
using std::vector;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void plotchart();
    void updateMessage(QString& message);
    void updateMessage(QString&& message);
    void updateChart(const double& x, const double& y);
private:
    Ui::MainWindow *ui;   
    QChart* m_chart = new QChart();
    QLineSeries* m_series = new QLineSeries;
    QValueAxis* m_axisX = new QValueAxis();
    QValueAxis* m_axisY = new QValueAxis();
    double X_MIN;
    double X_MAX;
    double Y_MIN;
    double Y_MAX;
};
#endif // MAINWINDOW_H
