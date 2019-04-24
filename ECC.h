// author: Hanjun Li <hanjunl@andrew.cmu.edu>

#ifndef ECC_H_
#define ECC_H_

#include <stdlib.h>

#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
// #include <libscapi/include/circuits/ArithmeticCircuit.hpp>
#include <libscapi/include/infra/Measurement.hpp>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include <libscapi/include/primitives/Mersenne.hpp>
#include "ProtocolTimer.h"
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/infra/Common.hpp>
// #include <libscapi/include/primitives/Prg.hpp>
// #include "HashEncrypt.h"
#include <emmintrin.h>
#include <thread>

#define flag_print false
#define flag_print_timings true
#define flag_print_output true


using namespace std;
using namespace std::chrono;



// Function 0: Reconstruct from ECC a polynomial.
// -- Players run locally.
// -- using Gao02 algorithm
template <class FieldType>
class ECC{

private:
  vector<FieldType> _g0;
  vector<FieldType> _alpha;
  vector<bool> _alphaMask;
  TemplateField<FieldType>*_field;
  
  void trimZeroes(vector<FieldType>& a);

  void applyMask(vector<FieldType>& fullVec);

  void addToPolynomial(vector<FieldType>& p1, // input
                       vector<FieldType>& p2);

  void multToPolynomial(vector<FieldType>& p1,
                        vector<FieldType>& p2);

  void scaleToPolynomial(FieldType c,
                         vector<FieldType>& p);

  void dividePolynomial(vector<FieldType>& p1,
                        vector<FieldType>& p2,
                        vector<FieldType>& q, // quotient
                        vector<FieldType>& r);

  void extendedEuclideanPartial(vector<FieldType>& p1, // input
                                vector<FieldType>& p2, // input
                                vector<FieldType>& a1, // output
                                vector<FieldType>& a2, // output
                                vector<FieldType>& r,  // output
                                int minDegree);
  void computeG0(vector<FieldType>& alpha, // input
                 vector<FieldType>& g0);   // output

public:
  ECC();
  ~ECC();

  // debugging
  void printPolynomial(vector<FieldType>& p);

  // independent version: using supplied x values
  FieldType evalPolynomial(FieldType x,
                           vector<FieldType>& polynomial);

  void interpolate(vector<FieldType>& x, // input 
                   vector<FieldType>& y, // input
                   vector<FieldType>& polynomial);

  bool reconstruct(vector<FieldType>& alpha, // input
                   vector<FieldType>& g0,    // input
                   vector<FieldType>& code,  // input
                   int degree,
                   vector<FieldType>& polynomial);

  
  // the class version: using a fixed alpha vector as x values
  void setAlpha(vector<FieldType>& alpha);
  void setField(TemplateField<FieldType>* field);
  void setPartyInactive(int partyId);
  
  
  
  void interpolate(vector<FieldType>& y,
                   vector<FieldType>& polynomial);

  bool reconstruct(vector<FieldType>& code, // input
                   int degree,
                   vector<FieldType>& polynomial);

};
  
// ---------------- implementations ----------------
template <class FieldType>
ECC<FieldType>::ECC(){
}

template <class FieldType>
ECC<FieldType>::~ECC(){
}


template<class FieldType> 
void ECC<FieldType>::
trimZeroes(vector<FieldType>& a) {
  int i = a.size();
  FieldType Zero = *(_field->GetZero());
  while(i > 0 && a[i-1] == Zero) { i--; }
  a.resize(i);
  return;
}

template <class FieldType>
void ECC<FieldType>::
addToPolynomial(vector<FieldType>& p1, // input
                vector<FieldType>& p2){ // both input and output
  // store result in p2
  int p1_size = p1.size();
  int p2_size = p2.size();
  if(p1_size > p2_size){
    p2.resize( p1_size, *(_field->GetZero()));
  }

  for(int i=0; i<p1_size; i++){
    p2[i] += p1[i];
  }
  trimZeroes(p2);
  return;
}

template <class FieldType>
void ECC<FieldType>::
multToPolynomial(vector<FieldType>& p1,
                 vector<FieldType>& p2){

  int p1_size = p1.size();
  int p2_size = p2.size();
  vector<FieldType> tempProduct(p1_size + p2_size - 1, *(_field->GetZero()));

  for(int i=0; i<p1_size; i++){
    if(p1[i] == *(_field->GetZero()) ){
      continue;
    }
    for(int j=0; j<p2_size; j++){
      tempProduct[i+j] += p1[i] * p2[j];
    }
  }
    
  trimZeroes(tempProduct);
  p2 = tempProduct;
  return;
}

template <class FieldType>
void ECC<FieldType>::
scaleToPolynomial(FieldType c,
                  vector<FieldType>& p){
  int p_deg = p.size();
  for(int i=0; i<p_deg; i++){
    p[i] *= c;
  }
  return;
}

template <class FieldType>
void ECC<FieldType>::
dividePolynomial(vector<FieldType>& p1,
                 vector<FieldType>& p2,
                 vector<FieldType>& q, // quotient
                 vector<FieldType>& r){ // remainder
  r = p1;
  int p1Size = p1.size();
  int p2Size = p2.size();
  if(p1Size < p2Size){
    q.resize(0);
    r = p1;
    return;
  }
  if (p2Size == 0) {
    cerr << "ECC: dividing by zero (polynomail)" << endl;
    abort();
  }

  int qSize = p1Size - p2Size +1;
  q.resize(qSize);

  for(int i = p1Size - p2Size; i>=0; i--){
    FieldType topCoeff = r[ p2Size + i -1] / p2[p2Size-1];

    q[i] = topCoeff;
    vector<FieldType> xi(i+1, *(_field->GetZero()) );
    xi[i] = xi[i] - topCoeff;
    vector<FieldType> p2Tmp = p2;
    
    // r -= topCoeff * p2 * x^i
    multToPolynomial(xi, p2Tmp);
    addToPolynomial(p2Tmp, r);
  }
  trimZeroes(q);
  trimZeroes(r);
  return;
}

template <class FieldType>
void ECC<FieldType>::
extendedEuclideanPartial(vector<FieldType>& p1, // input
                         vector<FieldType>& p2, // input
                         vector<FieldType>& a1, // output
                         vector<FieldType>& a2, // output
                         vector<FieldType>& r,  // output
                         int minDegree){        // input
  // p1 is previous remainder
  // p2 is previous divider
  // assuming trimed
  int rDeg = p1.size()-1;
  if(rDeg < minDegree){
    // keep reminder == (a1 = 1, a2 = 0)
    a1.resize(1);
    a1[0] = *(_field->GetOne() );
    a2.resize(0);
    r = p1;
    return;
  }

  vector<FieldType> remainder, nextA1, nextA2;
  dividePolynomial(p2, p1, a1, remainder);
  extendedEuclideanPartial(remainder, p1, nextA1, nextA2, r, minDegree);

  // a1 = nextA2 - quotient * nextA1
  scaleToPolynomial( *(_field->GetZero())  - *(_field->GetOne()), a1);
  multToPolynomial(nextA1, a1);
  addToPolynomial(nextA2, a1); 
  a2 = nextA1;
  return;
}



// template <class FieldType>
// void extendedEuclidean(vector<FieldType>& p1,
//                        vector<FieldType>& p2,
//                        FieldType& a1,
//                        FieldType& a2){}



template <class FieldType>
FieldType ECC<FieldType>::
evalPolynomial(FieldType x, 
               vector<FieldType>& polynomial){
  int degree = polynomial.size() - 1;
  // assert(degree >= 0);
  FieldType result = polynomial[0];
  FieldType x_value = *(_field->GetOne());
  for(int i=0; i<degree; i++){
    x_value *= x;
    result += polynomial[i+1] * x_value;
  }
  
  return result;
}


template <class FieldType>
void ECC<FieldType>::
interpolate(vector<FieldType>& x, // input 
            vector<FieldType>& y, // input
            vector<FieldType>& polynomial){

  // assert(y.size() == x.size());
  int nPoints = y.size();
  int degree = nPoints - 1;
  vector<FieldType> result(nPoints,  *(_field->GetZero()));

  // ---- O(n^2) to compute all numerators ----
  vector< vector<FieldType> > numerator_before_i(nPoints);
  vector< vector<FieldType> > numerator_skip_i(nPoints);
  // fill-in numerator_before_i from left to right
  numerator_before_i[0] = vector<FieldType>(1, *(_field->GetOne()));
  for(int i=1; i<nPoints; i++){
    numerator_before_i[i] = vector<FieldType>( 2, *(_field->GetOne()) );
    numerator_before_i[i][0] =  *(_field->GetZero()) - x[i-1];
    multToPolynomial( numerator_before_i[i-1], numerator_before_i[i] );
  }
  // fill-in numerator_after_i from right to left
  numerator_skip_i[nPoints-1] = vector<FieldType>(1, *(_field->GetOne()));
  for(int i=nPoints-2; i>=0; i--){
    numerator_skip_i[i] = vector<FieldType>(2, *(_field->GetOne()));
    numerator_skip_i[i][0] =  *(_field->GetZero()) - x[i+1];
    multToPolynomial( numerator_skip_i[i+1], numerator_skip_i[i] );
  }
  // multiply before_i and after_i to get skip_i
  for(int i=0; i<nPoints; i++){
    multToPolynomial(numerator_before_i[i], numerator_skip_i[i]);
  }
  
  // ---- O(n^2) to compute all factors ----
  vector<FieldType> factor(nPoints);
  for(int i=0; i<nPoints; i++){
    FieldType denom = evalPolynomial( x[i], numerator_skip_i[i] );
    factor[i] = y[i] / denom;
  }
  
  // ---- O(n^2) add to get result ----
  for(int i=0; i<nPoints; i++){
    scaleToPolynomial( factor[i], numerator_skip_i[i] );
    addToPolynomial(numerator_skip_i[i], result);
  }

  polynomial = result;
  return;
}

template <class FieldType>
void ECC<FieldType>::
applyMask(vector<FieldType> &fullVec) {
  int nActive = 0;
  int nTotal = _alpha.size();
  for (int i = 0; i < nTotal; i++) {
    if (_alphaMask[i]) {
      fullVec[nActive++] = fullVec[i];
    }
  }
  fullVec.resize(nActive);
  return;
}

template <class FieldType>
void ECC<FieldType>::
interpolate(vector<FieldType>& y,
            vector<FieldType>& polynomial){
  // use _alpha as x.
  vector<FieldType> activeAlpha = _alpha;
  vector<FieldType> activeY = y;
  applyMask(activeAlpha);
  applyMask(activeY);
  
  interpolate(activeAlpha, activeY, polynomial);
  return;
}


template<class FieldType> 
void ECC<FieldType>::
setAlpha(vector<FieldType>& alpha){

  trimZeroes(alpha);
  _alpha = alpha;
  int nPoints = _alpha.size();
  _alphaMask.clear();
  _alphaMask.resize(nPoints, true);

  if(nPoints == 0){
    return;
  }

  computeG0(alpha, _g0);
  return;
}

template<class FieldType> 
void ECC<FieldType>::
setField(TemplateField<FieldType>* field){
  _field = field;
  return;
}

template <class FieldType>
void ECC<FieldType>::
computeG0(vector<FieldType> &alpha,
          vector<FieldType> &g0) {

  int nPoints = alpha.size();
  // assert(nPoints > 0)
  
  g0.clear();
  g0.resize(2, *(_field->GetOne()));
  g0[0] =  *(_field->GetZero()) - alpha[0];
  vector<FieldType> factor(2, *(_field->GetOne()));

  for(int i=1; i<nPoints; i++){
    factor[0] =  *(_field->GetZero()) - alpha[i];
    multToPolynomial(factor, _g0);
  }
  return;
}

template <class FieldType>
void ECC<FieldType>::
setPartyInactive(int partyId) {
  if (_alphaMask[partyId] == false) {
    // nothing to do.
    return;
  }
  // otherwise, set mask
  _alphaMask[partyId] = false;

  // and also recompute _g0
  vector<FieldType> newAlpha = _alpha;
  applyMask(newAlpha);
  if (newAlpha.size() > 0) {
    computeG0(newAlpha, _g0);
  } else {
    _g0.clear();
  }

  return;
}

template <class FieldType>
void ECC<FieldType>::
printPolynomial(vector<FieldType>& p){
  int deg = 0;
  for(auto e : p){
    cout << e << "*x^" << (deg++) << " ";
  }
  cout << endl;
  return;
}

// from Gao02
// -- for a fixed alpha, g0 = prod_i(x - alpha[i])
// -- build g1 = interpolate(alpha, code) w/ degree = n-1
// -- use extended euclidean to find u*g0 + v*g1 = g
//    -- stop when g has degree smaller than (n+k)/2
// -- corrected message is g / v
// -- fail if remainder is not 0.
template <class FieldType>
bool ECC<FieldType>::
reconstruct(vector<FieldType>& alpha, // input x of size n
            vector<FieldType>& g0,    // input g0, of size n+1
            vector<FieldType>& code,  // input y of size n
            int degree,               // T-1
            vector<FieldType>& polynomial){

  int nMessages = degree+1; // k
  int nCodes = code.size(); // n
  vector<FieldType> g1;
  interpolate(code, g1);

  // minDeg = ceil( (n + k)/2 )
  int minDeg = (nCodes + nMessages +1)/2;

  vector<FieldType> u, v, g, r;
  // find u*g0 + v*g1 = g, with deg(g) < minDeg
  extendedEuclideanPartial(g0, g1, u, v, g, minDeg);

  // find f*v + r = g
  dividePolynomial(g, v, polynomial, r);

  if(r.size() > 0){
    return false;
  }
  return true;
}

template <class FieldType>
bool ECC<FieldType>::
reconstruct(vector<FieldType>& code, // input
            int degree,
            vector<FieldType>& polynomial){

  vector<FieldType> activeAlpha = _alpha;
  vector<FieldType> activeCode = code;
  applyMask(activeAlpha);
  applyMask(activeCode);

  return reconstruct(activeAlpha, _g0, activeCode, degree, polynomial);
}


#endif /* BAPARTY_H_ */
