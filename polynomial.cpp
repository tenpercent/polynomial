#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using std::vector;
using std::endl;
using std::cin;
using std::cout;
using std::ostream;

double const min_double = 1e-8;

template <typename T>
class Polynomial {
  vector<T> coefficients;
public:
  explicit Polynomial (T const & zeroth_coefficient) {
    coefficients.push_back(zeroth_coefficient);
    coefficients.push_back(0);
  }

  Polynomial() {
    coefficients.push_back(0);
  }

  template<class Iterator>
  Polynomial (Iterator const begin, Iterator const end) {
    for (typename vector<T>::iterator it = begin; it != end; ++it) {
      coefficients.push_back (*it);
    }
    coefficients.push_back(0);
  };

  bool operator== (Polynomial const & second) const;
  bool operator!= (Polynomial const & second) const;
  Polynomial operator+ (Polynomial const & second) const;
  Polynomial & operator+= (Polynomial const & second);
  Polynomial operator- (Polynomial const & second) const;
  Polynomial & operator-= (Polynomial const & second);
  Polynomial operator* (Polynomial const & second) const;
  Polynomial & operator*= (Polynomial const & second);
  Polynomial operator* (T const & second) const;
  Polynomial & operator*= (T const & second);

  int Degree () const;
  T & operator[] (size_t const degree);
  const T & operator[] (size_t const degree) const;
  T operator() (double const point) const;
  Polynomial operator/ (Polynomial const & second) const;
  Polynomial operator% (Polynomial const & second) const;
  Polynomial operator, (Polynomial const & second) const;

  void make_proper ();

  class iterator;
  Polynomial<T>::iterator begin();
  Polynomial<T>::iterator end();
};

template<typename T>
class Polynomial<T>::iterator : 
  public std::iterator<std::input_iterator_tag, T> {
    typename vector<T>::iterator pointed_value;
  public:
    explicit iterator (typename vector<T>::iterator initialize_from);
    iterator& operator++ ();
    iterator operator++ (int);
    T& operator*();
    bool operator== (iterator const & second);
    bool operator!= (iterator const & second);
  };

template<typename T>
ostream& operator<< (std::ostream & second, Polynomial<T> const & poly);

template<typename T>
void initialize_polynomial_ver1 (Polynomial<T> & polynomial) {
  polynomial[0] = 1;
  polynomial[1] = 4;
}

template<typename T>
void initialize_polynomial_ver2 (Polynomial<T> & polynomial) {
  polynomial[0] = 1;
  polynomial[1] = 2;
  polynomial[2] = 1;
}

int main () {
  Polynomial<double> first_polynomial;
  Polynomial<double> second_polynomial;
  Polynomial<double> third_polynomial;

  initialize_polynomial_ver1(first_polynomial);
  initialize_polynomial_ver2(second_polynomial);
  cout << "first_polynomial: " << first_polynomial << endl;
  cout << "second_polynomial: " << second_polynomial << endl;

  cout << "test degree" << endl;
  cout << first_polynomial.Degree() << endl;
  cout << second_polynomial.Degree() << endl;

  cout << "test operator+" << endl;
  third_polynomial = first_polynomial + second_polynomial;
  cout << third_polynomial << endl;
  
  cout << "test operator-" << endl;
  third_polynomial = first_polynomial - second_polynomial;
  cout << third_polynomial << endl;
  
  cout << "test operator*" << endl;
  third_polynomial = first_polynomial * second_polynomial;
  cout << third_polynomial << endl;

  cout << "test operator/" << endl;
  third_polynomial = first_polynomial / second_polynomial;
  cout << third_polynomial << endl;

  cout << "test operator%" << endl;
  third_polynomial = first_polynomial % second_polynomial;
  cout << third_polynomial << endl;

  cout << "test call operator: first(3) = ..." << endl;
  cout << first_polynomial(3) << endl;

  cout << "test gcd" << endl;
  third_polynomial = (first_polynomial, second_polynomial);
  cout << third_polynomial << endl;
  third_polynomial = (second_polynomial, first_polynomial);
  cout << third_polynomial << endl;

  cout << "test operator[] on get" << endl;
  for (size_t i = 0; i < 10; ++i) {
    cout << "first[" << i << "] = " << first_polynomial[i] << endl;
  }

  cout << "test operator[] on set" << endl;
  first_polynomial[4] = 2;
  cout << first_polynomial << endl;
  initialize_polynomial_ver1 (first_polynomial);

  cout << "test constructor by coefficient" << endl;
  Polynomial<double> fourth_polynomial = Polynomial<double>(4.2);
  cout << fourth_polynomial << endl;

  cout << "test constructor by iterators" << endl;
  vector<double> test_vector;
  for (size_t i = 0; i < 4; ++i) {
    test_vector.push_back((rand() % 10) * .1);
  }
  fourth_polynomial = Polynomial<double>(test_vector.begin(), test_vector.end());
  cout << fourth_polynomial << endl;

  cout << "test operator+=" << endl;
  initialize_polynomial_ver1 (first_polynomial);
  first_polynomial += first_polynomial;
  cout << first_polynomial << endl;

  cout << "test operator-=" << endl;
  initialize_polynomial_ver1 (first_polynomial);
  first_polynomial -= first_polynomial;
  cout << first_polynomial << endl;

  cout << "test operator*=" << endl;
  initialize_polynomial_ver1(first_polynomial);
  first_polynomial *= first_polynomial;
  cout << first_polynomial << endl;

  cout << "test iterators" << endl;
  for (
    Polynomial<double>::iterator it = first_polynomial.begin();
    it != first_polynomial.end();
    ++it)
  {
    cout << *it << endl;
  }
  return 0;
}

template<typename T>
void Polynomial<T>::make_proper() {
  for (int i = coefficients.size() - 1; i >= 0; --i) {
    if (fabs(coefficients[i]) < min_double) {
      coefficients.pop_back();
    } else {
      break;
    }
  }
  coefficients.push_back(0);
}

template<typename T>
  int Polynomial<T>::Degree () const {
    int coef_iterator;
    for (coef_iterator = coefficients.size() - 2; coef_iterator >= 0; --coef_iterator) {
      if (coefficients[coef_iterator] != 0) {
        break;
      }
    }
    return coef_iterator;
  }

template<typename T>
  bool Polynomial<T>::operator== (Polynomial<T> const & second) const {
    if ((this->Degree()) != (second.Degree())) {
      return false;
    }
    return (
      std::equal(
        this->coefficients.begin(), 
        this->coefficients.end(), 
        second.coefficients.begin()
      )
    );
  }

template<typename T>
  bool Polynomial<T>::operator!= (Polynomial const & second) const {
    return !(*this == second);
  }

template<typename T>
  const T & Polynomial<T>::operator[] (size_t const degree) const {
    if (degree >= coefficients.size()) {
      return coefficients.back();
    }
    return coefficients[degree];
  }

template<typename T>
  T & Polynomial<T>::operator[] (size_t const degree) {
    while (coefficients.size() < degree + 2) {
      coefficients.push_back(0);
    }
    return coefficients[degree];
    // here we've got a problem:
    // if we set leading coefficient to 0
    // it does affect some properties like degree
    // but we can't control it immediately
  }

template<typename T>
  Polynomial<T>& Polynomial<T>::operator+= (Polynomial const & second) {
    for (size_t i = 0; i < coefficients.size(); ++i) {
      coefficients[i] += second[i];
    }
    if ((coefficients.size()) < (second.coefficients.size())) {
      for (size_t i = coefficients.size(); i < second.coefficients.size(); ++i) {
        coefficients.push_back(second[i]);
      }
    }
    make_proper();
    return *this;
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator+ (Polynomial const & second) const {
    Polynomial<T> answer (*this);
    answer += second;
    return answer;
  }

template<typename T>
  ostream& operator<< (std::ostream & second, Polynomial<T> const & poly) {
    for (int i = poly.Degree(); i > 1; --i) {
      if (poly[i] == 1) {
         second << "x^" << i << " + ";
      } else {
         second << poly[i] << "x^" << i << " + ";
      }
    }
    if (poly.Degree() >= 1) {
      if (poly[1] == 1) {
        second << "x + ";
      } else {
        second << poly[1] << "x + ";
      }
    }
    second << poly[0];
    return second;
  }

template<typename T>
  Polynomial<T> & Polynomial<T>::operator*= (T const & multiplicator) {
    for (size_t i = 0; i < coefficients.size(); ++i) {
      coefficients[i] *= multiplicator;
    }
    return *this;
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator* (T const & multiplicator) const {
    Polynomial<T> product (*this);
    product *= multiplicator;
    return product;
  }

template<typename T>
  Polynomial<T> & Polynomial<T>::operator*= (Polynomial const & multiplicator) {
    Polynomial<T> this_copy_used_for_multiplication (*this);
    Polynomial<T> this_copy_used_for_storage (*this);
    Polynomial<T> temporal_product;
    this_copy_used_for_storage *= multiplicator[0];
    for (int i = 0; i < multiplicator.coefficients.size(); ++i) {
      for (size_t j = this_copy_used_for_multiplication.coefficients.size() + i; j > i; --j) {
        this_copy_used_for_multiplication[j] = this_copy_used_for_multiplication[j - 1];
      }     
      this_copy_used_for_multiplication[i] = 0;
      temporal_product = this_copy_used_for_multiplication * multiplicator[i + 1];
      this_copy_used_for_storage += temporal_product;
    }
    *this = this_copy_used_for_storage;
    make_proper();
    return *this;
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator* (Polynomial const & multiplicator) const {
    Polynomial<T> this_copy (*this);
    this_copy *= multiplicator;
    return this_copy;
  }

template<typename T>
  T Polynomial<T>::operator() (double const point) const {
    T answer (0);
    T temporal_multiplicator (1);
    for (int i = 0; i <= Degree(); ++i) {
      answer += (temporal_multiplicator * coefficients[i]);
      temporal_multiplicator *= point;
    }
    return answer;
  }

template<typename T>
  Polynomial<T>& Polynomial<T>::operator-= (Polynomial const & second) {
    Polynomial<T> subtracted = second * (-1);

    *this += subtracted;
    return *this;
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator- (Polynomial const & second) const {
    Polynomial<T> answer (*this);
    answer -= second;
    return answer;
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator/ (Polynomial const & second) const {
    Polynomial<T> divisor = second;
    divisor.make_proper();

    if (divisor.Degree() == -1) {
      return *this;
    }
    Polynomial<T> divided (*this);
    Polynomial<T> result;
    Polynomial<T> current_step_result;

    T result_coef;
    int result_degree = 0;
    int divided_degree = divided.Degree();
    int const divisor_degree = divisor.Degree();

    while (divided_degree >= divisor_degree) {
      current_step_result = Polynomial(0);
      result_coef = divided[divided_degree] / divisor[divisor_degree];
      result_degree = divided_degree - divisor_degree;
      current_step_result[result_degree] = result_coef;
      divided -= current_step_result * divisor;
      result += current_step_result;
      divided_degree = divided.Degree();
    }
    return result;
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator% (Polynomial const & second) const {
    Polynomial<T> result_poly = *this - (*this / second) * second;
    return result_poly;
/*
    Polynomial<T> divisor = second;
    divisor.make_proper();

    if (divisor.Degree() == -1) {
      return *this;
    }
    Polynomial<T> divided (*this);
    Polynomial<T> result;
    Polynomial<T> current_step_result;

    T result_coef;
    int result_degree = 0;
    int divided_degree = divided.Degree();
    int const divisor_degree = divisor.Degree();

    while (divided_degree >= divisor_degree) {
      current_step_result = Polynomial(0);
      result_coef = divided[divided_degree] / divisor[divisor_degree];
      result_degree = divided_degree - divisor_degree;
      current_step_result[result_degree] = result_coef;
      divided -= current_step_result * divisor;
      result += current_step_result;
      divided_degree = divided.Degree();
    }
    return divided;
*/
  }

template<typename T>
  Polynomial<T> Polynomial<T>::operator, (Polynomial const & second) const {
    Polynomial<T> smaller_poly = (this->Degree() < second.Degree()) ? *this : second;
    Polynomial<T> bigger_poly = (this->Degree() >= second.Degree()) ? *this : second;
    while (smaller_poly.Degree() >= 0) {
      bigger_poly = bigger_poly % smaller_poly;
      std::swap(bigger_poly, smaller_poly);
    }

    bigger_poly.make_proper();
    T leading_coef = bigger_poly[bigger_poly.Degree()];
    bigger_poly *= 1. / leading_coef;
    
    return bigger_poly;
  }

// iterator stuff

template<typename T>
  typename Polynomial<T>::iterator Polynomial<T>::begin() {
    return Polynomial<T>::iterator (coefficients.begin());
  }

template<typename T>
  typename Polynomial<T>::iterator Polynomial<T>::end() {
    return Polynomial<T>::iterator (coefficients.end() - 1);
  }

template<typename T>
  Polynomial<T>::iterator::iterator (typename vector<T>::iterator initialize_from) :
    pointed_value (initialize_from)
  {}

template<typename T>
  typename Polynomial<T>::iterator & Polynomial<T>::iterator::operator++() {
    ++pointed_value;
    return *this;
  }

template<typename T>
  typename Polynomial<T>::iterator Polynomial<T>::iterator::operator++(int) {
    Polynomial<T>::iterator this_copy = *this;
    operator++();
    return this_copy;
  }

template<typename T>
  T & Polynomial<T>::iterator::operator* () {
    return *pointed_value;
  }

template<typename T>
  bool Polynomial<T>::iterator::operator== (typename Polynomial<T>::iterator const & second) {
    return this->pointed_value == second.pointed_value;
  }

template<typename T>
  bool Polynomial<T>::iterator::operator!= (typename Polynomial<T>::iterator const & second) {
    return (!(operator==(second)));
  }
