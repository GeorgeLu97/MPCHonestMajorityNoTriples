#ifndef BAPARTY_H_
#define BAPARTY_H_

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


// Questions:
// 1. Is ``run in parallel'' different from ``run in batch''?
//    -- for now, run in sequence.
// 2. Error checking during protocol

// TODO:
// 1. look into Vandermonde matrix for ECC (VDM class?)
// 2. read consensus and broadcast papers
// 3. Error Correction (reconstrcut)
// 4. fault_localization. (store messages/ simulate party interactions)

// Function 0: Reconstruct from ECC a polynomial.
// -- Players run locally.
template <class FieldType>
bool reconstruct(vector<FieldType>& alpha, // input
                 vector<FieldType>& code, // input
                 int degree,
                 vector<FieldType>& polynomial);

template <class FieldType>
FieldType evalPolynomial(FieldType x,
                         vector<FieldType>& polynomial);

// Functionality 1: concensus protocol by BGP92
// Functionality 2: broadcast protocol by CW92
// Functionality 3: Robust batched broadcast in Appendex of BTH
// -- note: roll my own fault_localization
// -- note: fault_localization usess both consensus and broadcast
template <class FieldType>
class BAParty : public Protocol, HonestMajority, MultiParty{
  
private:

  // -------- protocol properties --------
  int _myId;
  int _smallT;
  vector<int> _parties;
  vector<int> _dealers;
  vector<FieldType> alpha; // <-- the commonly known field element (id) for each party
  HIM<FieldType> __peMatrix;

  // -------- private functionalities (sub-protocols) --------
  // from Appendex A. in BTH paper:
  // PE_Broadcast (each player spreads 1 value)
  // -- Every party P_i sends x_i to every P_j
  // -- Every party P_i computes (x''_1, ..., x''_n) = M(x'_1, ..., x'_n)
  // -- Every party P_i sends x''_i to evey P_j
  // -- Every party P_i checks if the received values are equal
  // -- Every party outputs received x_1, ..., x_n
  // returns happy bit.
  bool peBroadcast(const FieldType elem, // input
                   vector<FieldType>& received_elems);

  // BroadcastForP (each player spreads l values)
  // - Repeat t times (k = 0, ..., t-1):
  // -- Every party P_i sets a happy_bit = 1.
  // -- Run PE_Broadcast ceil(l/t) times, at a time
  // -- Every P_i sends to every P_j its happy_bit.
  // -- Every P_i ``AND'' all received happy bits.
  // -- Run concensus() to decide if there are unhappy players.
  // -- Run fault_localization()
  // -- Remove faulty players and repeat.
  // robust: no return values.
  void broadcastForAll(const vector<FieldType>& input_elems, // input
                       vector< vector<FieldType> >& output_elems);

public:
  // -------- public functionalities --------
  // constructor and destructors
  BAParty();
  ~BAParty();
  
  // set protocol properties: (should be known before running!)
  void setParties(const vector<int>& participants);
  void setDealers(const vector<int>& dealers);
  void remainingParties(vector<int>& participants);

  // from BGP92: consensus() only on a bit (happiness).
  // -- TODO, fill in design
  bool consensus(bool b);
  
  // from CW92:
  // -- TODO, fill in design
  void broadcast(vector<FieldType>& msg, bool isSender);

  // Broadcast (each dealer (k intotal) spreads T values)
  // -- Every dealer expand T values into n by interpolating
  // -- Every dealer distribute 1 value to each party
  // -- Run BroadcastForP() for k values (1 from each dealer)
  // -- Every player reconstruct T values for each dealer w/ ECC
  void robustBatchBroadcast(vector<FieldType>& elems, bool isDealer);
};



// ---------------- implementations ----------------
template<class FieldType> 
void trimZeroes(vector<FieldType>& a) {
  int i = a.size();
  FieldType Zero = FieldType(0);
  while(a[i-1] == Zero && i > 0) { i--; }
  a.resize(i);
  return;
}

template <class FieldType>
void addToPolynomial(vector<FieldType>& p1, // input
                     vector<FieldType>& p2){ // both input and output
  // store result in p2
  int p1_size = p1.size();
  int p2_size = p2.size();
  if(p1_size > p2_size){
    p2.resize( p1_size, FieldType(0) );
  }

  for(int i=0; i<p1_size; i++){
    p2[i] += p1[i];
  }
  trimZeroes(p2);
  return;
}

template <class FieldType>
void multToPolynomial(vector<FieldType>& p1,
                      vector<FieldType>& p2){

  int p1_size = p1.size();
  int p2_size = p2.size();
  vector<FieldType> tempProduct(p1_size + p2_size - 1, FieldType(0));
  for(int i=0; i<p1_size; i++){
    for(int j=0; j<p2_size; j++){
      tempProduct[i+j] += p1[i] * p2[j];
    }
  }
    
  trimZeroes(tempProduct);
  p2 = tempProduct;
  return;
}

template <class FieldType>
void scaleToPolynomial(FieldType c,
                       vector<FieldType>& p){
  int p_deg = p.size();
  for(int i=0; i<p_deg; i++){
    p[i] *= c;
  }
  return;
}

template <class FieldType>
FieldType evalPolynomial(FieldType x, 
                         vector<FieldType>& polynomial){
  int degree = polynomial.size() - 1;
  // assert(degree >= 0);
  FieldType result = polynomial[0];
  FieldType x_value = FieldType(1);
  for(int i=0; i<degree; i++){
    x_value *= x;
    result += polynomial[i+1] * x_value;
  }
  
  return result;
}

template <class FieldType>
void interpolate(vector<FieldType>& x, // input 
                 vector<FieldType>& y, // input
                 vector<FieldType>& polynomial){

  // assert(y.size() == x.size());
  int nPoints = y.size();
  int degree = nPoints - 1;
  vector<FieldType> result(nPoints, FieldType(0));

  // ---- O(n^2) to compute all numerators ----
  vector< vector<FieldType> > numerator_before_i(nPoints);
  vector< vector<FieldType> > numerator_skip_i(nPoints);
  // fill-in numerator_before_i from left to right
  numerator_before_i[0] = vector<FieldType>(1, FieldType(1));
  for(int i=1; i<nPoints; i++){
    numerator_before_i[i] = vector<FieldType>( 2, FieldType(1) );
    numerator_before_i[i][0] = FieldType(0) - x[i-1];
    multToPolynomial( numerator_before_i[i-1], numerator_before_i[i] );
  }
  // fill-in numerator_after_i from right to left
  numerator_skip_i[nPoints-1] = vector<FieldType>(1, FieldType(1));
  for(int i=nPoints-2; i>=0; i--){
    numerator_skip_i[i] = vector<FieldType>(2, FieldType(1));
    numerator_skip_i[i][0] = FieldType(0) - x[i+1];
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
  
  // add to get result
  for(int i=0; i<nPoints; i++){
    scaleToPolynomial( factor[i], numerator_skip_i[i] );
    addToPolynomial(numerator_skip_i[i], result);
  }

  polynomial = result;
  return;
}



// TODO: do real error correction.
// -- for now, a stub that interpolate the first T points.
template <class FieldType>
bool reconstruct(vector<FieldType>& alpha, // input x of size n
                 vector<FieldType>& code, // input y of size n
                 int degree, // T-1
                 vector<FieldType>& polynomial){
  int nPoints = degree+1;
  vector<FieldType> x(alpha.begin(), alpha.begin()+nPoints);
  vector<FieldType> y(code.begin(), code.begin()+nPoints);
  
  // TODO: change this
  interpolate(x, y, polynomial);

  return true;
}

// -------- private functions --------
template <class FieldType>
bool BAParty<FieldType>::
peBroadcast(const FieldType elem, // input
            vector<FieldType>& received_elems){
  return true;
  
}

template <class FieldType>
void BAParty<FieldType>::
broadcastForAll(const vector<FieldType>& input_elems, // input
                vector< vector<FieldType> >& output_elems){
  return;
}

// -------- constructor & destructor --------
// constructor and destructors
template <class FieldType>
BAParty<FieldType>::BAParty(){
}


template <class FieldType>
BAParty<FieldType>::~BAParty(){
}


// -------- public functionalities --------
// set protocol properties: (should be known before running!)
template <class FieldType>
void BAParty<FieldType>::
setParties(const vector<int>& participants){
  return;
}

template <class FieldType>
void BAParty<FieldType>::
setDealers(const vector<int>& dealers){
  return;
}

template <class FieldType>
void BAParty<FieldType>::
remainingParties(vector<int>& participants){
  return;
}

// from BGP92: consensus() only on a bit (happiness).
// -- TODO, fill in design
template <class FieldType>
bool BAParty<FieldType>::
consensus(bool b){
  return true;
}
  
// from CW92:
// -- TODO, fill in design
template <class FieldType>
void BAParty<FieldType>::
broadcast(vector<FieldType>& msg, bool isSender){
  return;
}

// Broadcast (each dealer (k intotal) spreads T values)
// -- Every dealer expand T values into n by interpolating
// -- Every dealer distribute 1 value to each party
// -- Run BroadcastForP() for k values (1 from each dealer)
// -- Every player reconstruct T values for each dealer w/ ECC
template <class FieldType>
void BAParty<FieldType>::
robustBatchBroadcast(vector<FieldType>& elems, bool isDealer){
  return;
}


#endif /* BAPARTY_H_ */
