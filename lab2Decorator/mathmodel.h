#ifndef MATHMODEL_H
#define MATHMODEL_H

#include "linearalgebra.h"
#include <math.h>

class IMathModel
{
protected:
    char name = 'A';
    long double sampIncrement, t_st, t_fin;
    matrix resMatrix;
    vector X0; // vector of initial condition
    int N;
public:
    IMathModel():sampIncrement(1e-1), t_st(0.0), t_fin(1.0),N(0){ resMatrix.resize(0, 0);}
    inline virtual vector getInitialConditions() const {return X0;}
    inline int getOrder()  {return X0.size();}
    inline virtual long double getSampIncrement() const {return this->sampIncrement;}
    inline virtual long double getT_st() const {return this->t_st;}
    inline virtual long double getT_fin() const {return this->t_fin;}
    inline virtual void setT_st(long double t_st) {this->t_st = t_st;}
    inline virtual void setT_fin(long double t_fin) {this->t_fin = t_fin;}
    inline virtual void setSampIncrement(long double sampIncrement) {this->sampIncrement = sampIncrement;}
    inline matrix getResult() const {return this->resMatrix;}
    virtual void getRP(const vector &X, long double t, vector &Y) const = 0; // X - current condition vector, Y - vector of increments
    char getName() const {return this->name;}
    int setX0byIndex(const int &index, const long double &value);
    inline long double getX0byIndex(const int &index) {return X0[index];}
    virtual void prepareResult();
    virtual void clearResult();
    virtual void preStart();
    virtual void addResult(const vector &addVector, const long double &t);
};

#endif // MATHMODEL_H
