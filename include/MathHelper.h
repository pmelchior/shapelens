#ifndef SHAPELENS_MATHHELPER_H
#define SHAPELENS_MATHHELPER_H

namespace shapelens {

  /// Templated integer powers: \f$x^n\f$.
  template <class T> 
  inline  T pow_int(const T& x, int n) {
    T x2,x3,x4;
    switch (n) {
    case 0: return T(1);
    case 1: return x;
    case 2: return x*x;
    case 3: return x*x*x;
    case 4: x2 = x*x; return x2*x2;
    case 5: x2 = x*x; return x2*x2*x;
    case 6: x3 = x*x*x; return x3*x3;
    case 7: x3 = x*x*x; return x3*x3*x;
    case 8: x2 = x*x; x4 = x2*x2; return x4*x4;
    case 9: x3 = x*x*x; return x3*x3*x3;
    case 10:x2 = x*x; x4 = x2*x2; return x4*x4*x2;
    default: x2 = T(1); // dummy variable for recursive calls
      while (n > 10) {
	x2 *= pow_int(x,10);
	n-=10;
      }
      return x2*pow_int(x,n); 
    }
  }
  /// Shorthands for \p pow_int(x,2) for small integer powers
  template <class T> inline T pow2(const T& x) { return pow_int(x,2); }
  /// Shorthands for \p pow_int(x,3) for small integer powers
  template <class T> inline T pow3(const T& x) { return pow_int(x,3); }
  /// Shorthands for \p pow_int(x,4) for small integer powers
  template <class T> inline T pow4(const T& x) { return pow_int(x,4); }
  /// Shorthands for \p pow_int(x,5) for small integer powers
  template <class T> inline T pow5(const T& x) { return pow_int(x,5); }
  /// Shorthands for \p pow_int(x,6) for small integer powers
  template <class T> inline T pow6(const T& x) { return pow_int(x,6); }
  /// Shorthands for \p pow_int(x,7) for small integer powers
  template <class T> inline T pow7(const T& x) { return pow_int(x,7); }
  /// Shorthands for \p pow_int(x,8) for small integer powers
  template <class T> inline T pow8(const T& x) { return pow_int(x,8); }
  /// Shorthands for \p pow_int(x,9) for small integer powers
  template <class T> inline T pow9(const T& x) { return pow_int(x,9); }
  /// Shorthands for \p pow_int(x,10) for small integer powers
  template <class T> inline T pow10(const T& x) { return pow_int(x,10); }

  /// Factorial of \p n.
  /// \b CAUTION: this overflows for \p n > 10!!!
  inline unsigned long factorial(int n) {
    unsigned long f = 1;
    for (int m=2; m <= n; m++)
      f *= m;
    return f;
  }

  /// Binomial coefficient \f$\begin{pmatrix}n\\m\end{pmatrix}\f$.
  /// \b CAUTION: This overflows for integer arguments larger 10!!!
  inline unsigned long binomial(int n, int m) {
    return factorial(n)/(factorial(m)*factorial(n-m));
  }
}


#endif
