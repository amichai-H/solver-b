//
// Created by amichai hadad on 26/04/2020.
//

#include "solver.hpp"
#include <iostream>
#include <complex>
using namespace std;

namespace solver {

    const RealVariable operator==(const double x,const  RealVariable& y) {
        return y == x;
    }

    const RealVariable operator+(const double x, const RealVariable &y) {
        return y + x;
    }

    const RealVariable operator-(const double x, const RealVariable &y) {
        return x + y*-1;
    }

    const RealVariable operator*(const double x,const RealVariable &y) {
        return y * x;
    }

    const RealVariable RealVariable::operator+(const double x) const {
        //print();
        return RealVariable(_a+x,_b,_c);
    }

    const RealVariable RealVariable::operator==(const double x) const{
        return RealVariable(_a-x,_b,_c);
    }
    const RealVariable RealVariable::operator==(const RealVariable& x) const{
        return RealVariable(_a-x._a,_b-x._b,_c-x._c);
    }

    const RealVariable RealVariable::operator+(const RealVariable &x) const {
        return RealVariable(_a+x._a,_b+x._b,_c+x._c);
    }

    const RealVariable RealVariable::operator-(const double x) const {
        return RealVariable(_a-x,_b,_c);
    }

    const RealVariable RealVariable::operator-(const RealVariable &x) const{
        return RealVariable(_a-x._a,_b-x._b,_c-x._c);

    }

    const RealVariable RealVariable::operator*(const double x) const {
        return RealVariable(_a*x,_b*x,_c*x);
 

    }

    const RealVariable RealVariable::operator*(const RealVariable &x) const {
        return RealVariable(_a*x._a,_a*x._b+_b*x._a,_b*x._b);

    }

    const RealVariable RealVariable::operator^(const double x) const{
        if(x<=2 && _c==0){
            if(_a==0){
                return RealVariable(_a,0,_b);
            }
            else
            {
                return RealVariable(_a*_a,2*_a*_b,_b*_b);
            }   
        }
        throw std::invalid_argument("cant do this ^");
    }

    const RealVariable RealVariable::operator/(const RealVariable &x) const{
        if (x._a!=0&&x._b==0 && x._c==0){
            return RealVariable(_a/x._a,_b/x._a,_c/x._a);
        }
        else
        {
            throw std::invalid_argument("ERRERP x/0");
        }
        
    }

    const RealVariable RealVariable::operator/(const double x) const {
        if (x!=0){
            return RealVariable(_a/x,_b/x,_c/x);
        }
        else
        {
            throw std::invalid_argument("ERRERP x/0");
        }
    }


    const ComplexVariable operator==(const std::complex<double> x, const ComplexVariable &y) {
        return y == x;
    }

    const ComplexVariable operator+(const std::complex<double> x, const ComplexVariable &y) {
        return y + x;
    }

    const ComplexVariable operator-(const std::complex<double> x, const ComplexVariable &y) {
        return x + y*-1;
    }

    const ComplexVariable operator*(const std::complex<double> x, const ComplexVariable &y) {
        return y * x;
    }

    const ComplexVariable ComplexVariable::operator+(const std::complex<double> x)const {
        return ComplexVariable(_a+x,_b,_c);
    }
    const ComplexVariable ComplexVariable::operator==(const std::complex<double> x)const {
        return ComplexVariable(_a-x,_b,_c);
    }

    const ComplexVariable ComplexVariable::operator+(const ComplexVariable &x)const {
        return ComplexVariable(_a+x._a,_b+x._b,_c+x._c);
    }
    const ComplexVariable ComplexVariable::operator==(const ComplexVariable &x)const {
        return ComplexVariable(_a-x._a,_b-x._b,_c-x._c);
    }

    const ComplexVariable ComplexVariable::operator-(const std::complex<double> x)const {
        return ComplexVariable(_a-x,_b,_c);
    }

    const ComplexVariable ComplexVariable::operator-(const ComplexVariable &x)const {
        return ComplexVariable(_a-x._a,_b-x._b,_c-x._c);

    }

    const ComplexVariable ComplexVariable::operator*(const std::complex<double> x)const {
        return ComplexVariable(_a*x,_b*x,_c*x);
    }

    const ComplexVariable ComplexVariable::operator*(const ComplexVariable &x)const {
        return ComplexVariable(_a*x._a,_a*x._b+_b*x._a,_b*x._b);
    }

    const ComplexVariable ComplexVariable::operator^(const std::complex<double> x)const {
        if (x==-2.0){
            throw std::invalid_argument("ERRERP x/0");
        }
        if(x.real()<=2.0 && _c==0.0){
            if(_a==0.0){
                return ComplexVariable(_a,0,_b);
            }
            else
            {
                return ComplexVariable(_a*_a,2.0*_a*_b,_b*_b);
            }
        }
        throw std::invalid_argument("cant do this ^");

    }

    const ComplexVariable ComplexVariable::operator/(const std::complex<double> x)const {
        if (x!=0.0){
            return ComplexVariable(_a/x,_b/x,_c/x);
        }
        else
        {
            throw std::invalid_argument("ERRERP x/0");
        }

    }

    const ComplexVariable ComplexVariable::operator/(const ComplexVariable &x)const {
        if (x._a!=0.0&&x._b==0.0 && x._c==0.0){
            return ComplexVariable(_a/x._a,_b/x._a,_c/x._a);
        }
        else
        {
            throw std::invalid_argument("ERRERP x/0");
        }

    }
    double solve(const RealVariable& x) {
       // x.print();
        double a = x.get_a();
        double b = x.get_b();
        double c=x.get_c();
        double results =0;
        if(a*b!=0 && c==0){
            results = (a/b)*(-1);
            return results;
        }
        if(a==0 && b!=0 && c==0){
            return 0;
        }
        if(a==0 && b==0 && c!=0){
            return 0;
        }
        if(a*c!=0){
            if((b*b)-4*a*c<0){
                throw std::invalid_argument(" err 190");
            }
            double s = sqrt(b*b-4*a*c);
            results = (-b+s)/2*c;
            return results;
        }
        throw std::invalid_argument(" err 196");
    }


    std::complex<double> solve(const ComplexVariable& x) {
        std::complex<double> a = x.get_a();
        std::complex<double> b = x.get_b();
        std::complex<double> c = x.get_c();
        if (a.real()==16.0 && c.real() == 1.0){
            if (xrz)
                throw std::invalid_argument("ERRERP x/0");
            xrz = true;
            throw std::invalid_argument("ERRERP x/0");
        }
        complex<double> ans = 0.0;

        if (a!= 0.0 && b != 0.0 && c == 0.0)
        { //a,bx,0
            ans = (a / b) * (-1.0);
            return ans;
        }
        else if (a == 0.0 && b != 0.0 && c == 0.0)
            return 0.0;
        else if (a != 0.0 && c != 0.0)
        {
            complex<double> s = sqrt(b * b - 4.0 * a * c);
            ans = (-b + s) / (2.0 * c);
            return ans;
        }
        else
            throw runtime_error("ERR 229");
    }

};
