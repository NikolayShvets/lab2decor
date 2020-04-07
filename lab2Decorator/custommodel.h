#ifndef CUSTOMMODEL_H
#define CUSTOMMODEL_H

#include "mathmodel.h"
#include "simplealgorithms.h"

class satellite : public IMathModel
{
private:
    const long double Re{6371000.0L};
    const long double nu{398600.436e9}; //грав параметр
public:
    satellite();
    void preStart() override;
    void getRP(const vector &X,long double t, vector &Y) const override;
};

class moon : public IMathModel
{
private:
    const long double Re{6371000.0L};
    const long double nu{398600.436e9}; //грав параметр
public:
    moon();
    void preStart() override;
    void getRP(const vector &X,long double t, vector &Y) const override;
};

class SomethingDecoratesSatellite : public IMathModel
{
protected:
    IMathModel *model;
public:
    SomethingDecoratesSatellite(IMathModel *model);
    vector getInitialConditions() const override {return model->getInitialConditions();}
    void setT_st(long double t_st) override {model->setT_st(t_st);}
    void setT_fin(long double t_fin) override {model->setT_fin(t_fin);}
    void setSampIncrement(long double sampIncrement) override {model->setSampIncrement(sampIncrement);}
    long double getT_st() const override {return model->getT_st();}
    long double getT_fin() const override {return model->getT_fin();}
    long double getSampIncrement() const override {return model->getSampIncrement();}
    void clearResult() override{model->clearResult();}
    void prepareResult() override {model->prepareResult();}
    void addResult(const vector &addVector, const long double &t) override {model->addResult(addVector, t);}
};

class MoonDecoratesSatellite : public SomethingDecoratesSatellite
{
private:
    const long double nu{4902.800e9};
    moon *m;
public:
    MoonDecoratesSatellite(IMathModel *model, moon *m);
    void getRP(const vector &X, long double t, vector &Y) const override;
};

#endif // CUSTOMMODEL_H
