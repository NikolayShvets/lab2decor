#include "custommodel.h"

satellite::satellite():IMathModel(){
    name = 'S';
    X0.resize(6);
    X0[0] = 0.9;  //big omega
    X0[1] = 1.2;  //naklonenie
    X0[2] = 0.85;  //pericenter latitude
    X0[3] = 19440000.0; //big axis
    X0[4] = 0.0; //e
    X0[5] = 0.8;  //true anomal

    X0[3] += Re;
    preStart();
}

void satellite::preStart() {
    this->X0 = simpleAlgorithms::orbitToXYZ(X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], nu);
}

void satellite::getRP(const vector &X,long double t, vector &Y) const{
    Y.resize(X.size() + 3);
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    long double ro = sqrt(pow(X[0],2.) + pow(X[1], 2.) + pow(X[2], 2.));
    Y[3] = -nu*X[0]/pow(ro,3.) + Y[6];
    Y[4] = -nu*X[1]/pow(ro,3.) + Y[7];
    Y[5] = -nu*X[2]/pow(ro,3.) + Y[8];
}

moon::moon():IMathModel(){
    name = 'M';
    X0.resize(6);
    X0[0] = 0.9;  //big omega
    X0[1] = 0.08726646;  //naklonenie
    X0[2] = 0.85;  //pericenter latitude
    X0[3] =  384748000.0; //big axis
    X0[4] = 0.0549; //e
    X0[5] = 0.8;  //true anomal

    X0[3] += Re;
    preStart();
}

void moon::preStart(){
    this->X0 = simpleAlgorithms::orbitToXYZ(X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], nu);
}

void moon::getRP(const vector &X,long double t, vector &Y) const{
    Y.resize(X.size());
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    long double ro = sqrt(pow(X[0],2.) + pow(X[1], 2.) + pow(X[2], 2.));
    Y[3] = -nu*X[0]/pow(ro,3.);
    Y[4] = -nu*X[1]/pow(ro,3.);
    Y[5] = -nu*X[2]/pow(ro,3.);
}

SomethingDecoratesSatellite::SomethingDecoratesSatellite(IMathModel *model):IMathModel(){
    this->model = model;
}

MoonDecoratesSatellite::MoonDecoratesSatellite(IMathModel *model, moon *m):SomethingDecoratesSatellite(model)
{
    this->m = m;
    X0.resize(6);
}

void MoonDecoratesSatellite::getRP(const vector &X, long double t, vector &Y) const
{
    long double minTime{999999999.0};
    int minTimeIndex;
    for(int i = 0; i < m->getResult().rowsCount() - 1; i++){
        if(std::abs(t - m->getResult()(i, m->getResult().colsCount() - 1)) < minTime){
            minTimeIndex = i;
            minTime = m->getResult()(i, m->getResult().colsCount() - 1);
        }
    }
    Y.resize(X.size() + 3);
    long double ro = sqrt(pow(X[0] - m->getResult()(0,minTimeIndex),2.) + pow(X[1] - m->getResult()(1,minTimeIndex),2.) + pow(X[2] - m->getResult()(2,minTimeIndex), 2.));
    Y[6] = nu*X[0]/pow(ro,3.);
    Y[7] = nu*X[1]/pow(ro,3.);
    Y[8] = nu*X[2]/pow(ro,3.);
    /*как теперь эту добавку учитывать?*/
    model->getRP(X, t, Y);
}
