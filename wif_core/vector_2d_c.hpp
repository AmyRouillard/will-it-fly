#ifndef __VECTOR_2D_C_INCLUDED__
#define __VECTOR_2D_C_INCLUDED__
#include <cmath>
#include <ostream>

namespace wif_core
{

class vector_2d_c
{
public:
    double x;
    double y;
    vector_2d_c(double const, double const); /// constuctior with x,y
    vector_2d_c();
    //vector_2d_c(vector_2d_c &); /// copy contuctor
    vector_2d_c operator + (vector_2d_c const) const; /// normal addition
    vector_2d_c operator -() const; /// returns negative vector
    vector_2d_c operator - (vector_2d_c const) const; /// normal substaction
    vector_2d_c operator * (double const) const;
    vector_2d_c operator / (double const) const;
    void operator += (vector_2d_c const); /// modifies vector, adds seconds vector
    void operator -= (vector_2d_c const); /// modifies vector, substacts seconds vector
    void operator *= (double const); /// modifies vector , multiplices
    void operator /= (double const); /// modifies vector , multiplices
    double get_length() const;/// lenght of vector
    double get_length_sq() const;/// same as pow(get_length(),2) but faster
    double get_angle() const; ///angle with x-axis
    double cross(vector_2d_c const) const; ///magnetude of 3d cross product. Used in intersection math.
}; // class


vector_2d_c operator * (double, vector_2d_c);
std::ostream & operator << (std::ostream & , const vector_2d_c & );


} // namespace wif_core

#endif // __VECTOR_2D_C_INCLUDED__
