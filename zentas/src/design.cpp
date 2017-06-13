#include <utility>
#include <iostream>

/* base clustering class */
class Base{
  public:
    int x;
    /* function specific to metric or data */
    virtual void f() = 0;
    Base(int x_):x(x_){}
};

/* optimised level 1*/
class OPT1 : public Base{
  public:
    int x;
    OPT1(int x_): Base(1), x(x_) {}
  
};

/* optimised level 2*/
class OPT2 : public Base{
  public:
    double y;
    double z;
    OPT2(double y_, double z_):Base(2), y(y_), z(z_) {}  
};

/* metric or data class */
class MET1{
  public:
    void print(){
      std::cout << "I am a MET1" << std::endl;
    }
};

/* here, tell Base (via TOPT) how TMET is used for f */
template <class TMET, class TOPT>
class X : public TOPT{
  public:
    
    template <class ... T>
    X(T ... t): TOPT(t...){}
    
    virtual void f() override final{
      TMET b;
      b.print();
    }
};


void test(){
  X<MET1, OPT1> p(13);
  X<MET1, OPT2> q(13,11.1);
  q.f();
}
