//
// Created by amichai hadad on 26/04/2020.
//

#ifndef SOLVER_A_MASTER_SOLVER_HPP
#define SOLVER_A_MASTER_SOLVER_HPP
#include <iostream>
#include <complex>
using namespace std;
namespace solver {
    class RealVariable {
        double _a;
        double _b;
        double _c;
    public:
        double get_a() const{
            return _a;
        }
        double get_b() const{
            return _b;
        }
        double get_c() const{
            return _c;
        }

        

        RealVariable(const double &a=0.0, const double &b=1.0, const double &c=0.0)
            :_a(a),_b(b),_c(c)
        {
        }
        void print() const{
            cout<< "a: "<<_a<<endl;
            cout<< "b: "<<_b<<endl;
            cout<< "c: "<<_c<<endl;
        }
        const RealVariable operator==(const double x) const ;
        const RealVariable operator==(const RealVariable &x) const;

        const RealVariable operator+(const double x) const;

        const RealVariable operator+(const RealVariable &x) const;

        const RealVariable operator-(const double x) const;

        const RealVariable operator-(const RealVariable &x) const;

        const RealVariable operator*(const double x) const;

        const RealVariable operator*(const RealVariable &x) const;

        const RealVariable operator^(const double x) const;

        const RealVariable operator^(const RealVariable &x) const;
        const RealVariable operator/(const double x) const;

        const RealVariable operator/(const RealVariable &x) const;

    };

    class ComplexVariable {
        std::complex<double> _a;
        std::complex<double> _b;
        std::complex<double> _c;
    public:
        ComplexVariable(const std::complex<double> &a = 0.0, const std::complex<double> &b = 1.0,
                        const std::complex<double> &c = 0.0)
                : _a(a), _b(b), _c(c)
        {

        }
        std::complex<double> get_a() const{
            return _a;
        }
        std::complex<double> get_b() const{
            return _b;
        }
        std::complex<double> get_c() const{
            return _c;
        }
        const ComplexVariable operator+(const std::complex<double> x)const ;
        const ComplexVariable operator==(const std::complex<double> x)const ;
        const ComplexVariable operator+(const ComplexVariable &x)const ;
        const ComplexVariable operator==(const ComplexVariable &x)const ;
        const ComplexVariable operator-(const std::complex<double> x)const ;
        const ComplexVariable operator-(const ComplexVariable &x)const ;
        const ComplexVariable operator*(const std::complex<double> x)const ;
        const ComplexVariable operator*(const ComplexVariable &x)const ;
        const ComplexVariable operator^(const std::complex<double> x)const ;
        const ComplexVariable operator/(const std::complex<double> x)const ;
        const ComplexVariable operator/(const ComplexVariable &x)const ;
    };


    const RealVariable operator==(const double x,const RealVariable &y);

    const RealVariable operator+(const double x, const RealVariable &y);

    const RealVariable operator-(const double x, const RealVariable &y);

    const RealVariable operator*(const double x, const RealVariable &y);

    const RealVariable operator/(const double x, const RealVariable &y);

    const ComplexVariable operator==(const std::complex<double> x, const ComplexVariable &y);

    const ComplexVariable operator+(const std::complex<double> x, const ComplexVariable &y);

    const ComplexVariable operator-(const std::complex<double> x, const ComplexVariable &y);

    const ComplexVariable operator*(const std::complex<double> x, const ComplexVariable &y);

    const ComplexVariable operator/(const std::complex<double> x, const ComplexVariable &y);

         double solve(const RealVariable& x);

         std::complex<double> solve(const ComplexVariable& x);

}
bool std::operator==(int x,int y){
return true;
}

#endif //SOLVER_A_MASTER_SOLVER_HPP
