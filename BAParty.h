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
FieldType evalPolynomial(vector<FieldType>& polynomial,
                         FieldType x);

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



template <class FieldType>
FieldType evalPolynomial(vector<FieldType>& polynomial,
                         FieldType x){
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
  // polynomial.resize(nPoints);
  
  lagrange(x, y, polynomial);
  
  return;
}



// TODO: do real error correction.
// -- for now, a stub that interpolate the first T points.
template <class FieldType>
bool reconstruct(const vector<FieldType>& alpha, // input x of size n
                 const vector<FieldType>& code, // input y of size n
                 int degree, // T-1
                 vector<FieldType>& polynomial){
  int nPoints = degree+1;
  vector<FieldType> x(alpha.begin(), alpha.begin()+nPoints);
  vector<FieldType> y(code.begin(), code.begin()+nPoints);
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
