#include <iostream>
#include "integrator.h"
#include "custommodel.h"
using namespace std;

int main()
{
    IMathModel *s;
    dormandPrinceIntgrator i;
    long double precision{}, t0{}, tk{}, dt{};
    bool isLunar{false};
    cout<<"input precision: "<<endl;
    cin>>precision;
    i.setPrecision(precision);
    cout<<"input t0, tk, dt: "<<endl;
    cin>>t0>>tk>>dt;
    cout<<"lunar?"<<endl;
    cin>>isLunar;
    s = new satellite;
    if(isLunar){
        moon *m = new moon;
        m->setT_st(t0);
        m->setT_fin(tk);
        m->setSampIncrement(dt);
        i.run(m);
        s = new MoonDecoratesSatellite(s,m);
    }
    s->setT_st(t0);
    s->setT_fin(tk);
    s->setSampIncrement(dt);
    i.run(s);
    return 0;
}
