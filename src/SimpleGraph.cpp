#include <iostream>
// Bugfix: Disable C++11 feature macros on Mac/libstdc++ to fool Qt into not using C++11 headers
#ifdef __GLIBCXX__
#ifdef __APPLE__
#undef __cplusplus
#define __cplusplus 200303L
#undef __GXX_EXPERIMENTAL_CXX0X__
#endif
#endif

#include <QtGui>
#include <QWidget>
#include <QApplication>
#include <algorithm>
#include <QCoreApplication>

#include "SimpleGraph.h"
#undef main

class MyWidget : public QWidget {
    public:
        static MyWidget & getInstance();
        void DrawGraph(SimpleGraph & userGraph);

    protected:
        SimpleGraph userGraph;
        void paintEvent(QPaintEvent *event);
};

MyWidget& MyWidget::getInstance() {
    static MyWidget instance;
    return instance;
}


void MyWidget::DrawGraph(SimpleGraph & userGraph) {
    this->userGraph = userGraph;
}

void MyWidget::paintEvent(QPaintEvent *event) {
    Q_UNUSED(event);
    static QTime last_run = QTime::currentTime();
    while (last_run.msecsTo(QTime::currentTime()) < 20);
    QPainter painter(this);
    if (!userGraph.nodes.empty()) {
        SimpleGraph copy = userGraph;
        auto getX = [](const Node& a, const Node& b) {return a.x < b.x;};
        double maxX = std::max_element(copy.nodes.begin(), copy.nodes.end(), getX)->x;
        double minX = std::min_element(copy.nodes.begin(), copy.nodes.end(), getX)->x;

        auto getY = [](const Node& a, const Node& b) { return a.y < b.y;};
        double maxY = std::max_element(copy.nodes.begin(), copy.nodes.end(), getY)->y;
        double minY = std::min_element(copy.nodes.begin(), copy.nodes.end(), getY)->y;

        auto scale = [maxX, maxY, minX, minY](const Node& a) -> Node {
            return {(a.x - minX) * 590 / (-minX + maxX) + 5, (a.y - minY) * 590 / (-minY + maxY) + 5};
        };
        painter.setBrush(Qt::blue);
        std::transform(copy.nodes.begin(), copy.nodes.end(), copy.nodes.begin(), scale);
        for (Edge & e : copy.edges) {
            Node start = copy.nodes[e.start];
            Node end = copy.nodes[e.end];
            painter.drawLine(start.x, start.y, end.x, end.y);
        }
        for (Node & n : copy.nodes) {
            if (n.type=='G'|| n.type=='C' ) painter.setBrush(Qt::red);
            //if (n.type=='A'|| n.type=='U' || n.type=='T' ) painter.setBrush(Qt::red);
            painter.drawEllipse(n.x - 10, n.y - 10, 20, 20);
        }
    }
    last_run = QTime::currentTime();
    update();
}

class WorkerThread : public QThread {
    void run() {
        int _userMain();
        _userMain();
    }
};

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    MyWidget & myWidget = MyWidget::getInstance();
    myWidget.resize(1200, 1200);
    myWidget.show();
    WorkerThread x;
    x.start();
    return app.exec();
}

void DrawGraph(SimpleGraph& userGraph) {
    MyWidget::getInstance().DrawGraph(userGraph);
}
