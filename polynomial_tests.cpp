template<typename T>
void initialize_polynomial_ver1(Polynomial<T> & polynomial) {
  polynomial[0] = 0;
  polynomial[1] = 1;
  polynomial[2] = 2;
}

template<typename T>
void initialize_polynomial_ver2(Polynomial<T> & polynomial) {
  polynomial[4] = 5;
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

  cout << "test operator[] on get" << endl;
  for (size_t i = 0; i < 10; ++i) {
    cout << "first[" << i << "] = " << first_polynomial[i];
  }

  cout << "test operator[] on set" << endl;
  first[4] = 2;
  cout << first_polynomial << endl;
  initialize_polynomial_ver1 (first_polynomial);

  cout << "test constructor by coefficient" << endl;
  Polynomial<double> fourth_polynomial = Polynomial(4.2);
  cout << fourth_polynomial << endl;

  cout << "test constructor by iterators" << endl;
  vector<double> test_vector;
  for (size_t i = 0; i < 4; ++i) {
    test_vector.push_back((rand() % 10) * .1);
  }
  fourth_polynomial = Polynomial(test_vector.begin(), test_vector.end());
  cout << fourth_polynomial << endl;

  cout << "test operator+=" << endl;
  initialize_polynomial_ver1(first_polynomial);
  first_polynomial += first_polynomial;
  cout << first_polynomial << endl;

  cout << "test operator-=" << endl;
  initialize_polynomial_ver1(first_polynomial);
  first_polynomial -= first_polynomial;
  cout << first_polynomial << endl;

  cout << "test operator*=" << endl;
  initialize_polynomial_ver1(first_polynomial);
  first_polynomial *= first_polynomial;
  cout << first_polynomial << endl;

  return 0;
}