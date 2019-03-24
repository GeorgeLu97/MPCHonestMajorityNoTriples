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
// 2. Error checking during protocol (e.g. timeout communication?)
// 3. How to agree on a constant HIM matrix?
// 4. send byte instead of bits?

// TODO:
// 2. read broadcast papers
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
class BAParty{
  
private:

  // -------- protocol properties --------
  int _myId;                    // my party id
  int _smallT;                  // number of assumed adversaries
  int _nThread;
  vector<bool> _dealersMask;    // mask to dealers
  vector<bool> _activeMask;     // mask to active players
  vector< shared_ptr<ProtocolPartyData> > _partySocket; // for communication
  vector<FieldType> _alpha; // <-- the commonly known elements to all
  vector<FieldType> _beta;  //     ..one for each party
  HIM<FieldType> _peMatrix; // common HIM matrix for PE_Broadcast
  
  // -------- private helper functions --------
  void exchangeBitWorker(bool sendBit, vector<bool>& recvBits,
                         int threadID, int nThreads);
  void exchangeBit(bool sendBit, vector<bool>& recvBits);

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
  // constructor and destructors
  BAParty();
  ~BAParty();
  
  // set protocol properties: (should be known before running!)
  // TODO: change depending on how outter protocol is implemented (ask George)
  void setParties(vector< shared_ptr<ProtocolPartyData> >& parties, int myId);
  void setHIM(HIM<FieldType>& common_HIM);
  void setAlphaBeta(vector<FieldType>& alpha, vector<FieldType>& beta);
  void setDealers(const vector<int>& dealers);
  void setSmallT(int smallT);
  void setNumThreads(int numThreads);
  void getRemainingParties(vector< shared_ptr<ProtocolPartyData> >& parties);


  // -------- public functionalities --------
  // from BGP92: consensus() only on a bit (happiness).
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
  
  // ---- O(n^2) add to get result ----
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
setParties(vector< shared_ptr<ProtocolPartyData> >& parties, int myId){
  int nParties = parties.size();
  _partySocket = parties;
  _myId = myId;
  // calculate defult number of advarseries
  _smallT = (nParties -1) / 3;
  // set default masks
  _dealersMask.clear();
  _dealersMask.resize(nParties, false);
  _activeMask.clear();
  _activeMask.resize(nParties, true);
  return;
}

template <class FieldType>
void BAParty<FieldType>::
setHIM(HIM<FieldType>& common_HIM){
  _peMatrix = common_HIM;
  return;
}

template <class FieldType>
void BAParty<FieldType>::
setAlphaBeta(vector<FieldType>& alpha, vector<FieldType>& beta){
  // assert(alpha.size() == _partySocket.size());
  // assert(beta.size() == _partySocket.size());
  _alpha = alpha;
  _beta = beta;
  return;
}

template <class FieldType>
void BAParty<FieldType>::
setDealers(const vector<int>& dealers){
  _dealersMask.clear();
  _dealersMask.resize( _partySocket.size(), false );
  for(auto d : dealers){
    // assert(d < _partySocket.size());
    _dealersMask[d] = true;
  }
  return;
}

template <class FieldType>
void BAParty<FieldType>::
setSmallT(int smallT){
  // assert( smallT < (_partySocket.size() -1) / 3 );
  _smallT = smallT;
  return;
}

template <class FieldType>
void BAParty<FieldType>::
setNumThreads(int numThreads){
  _nThread = numThreads;
  return;
}

template <class FieldType>
void BAParty<FieldType>::
getRemainingParties(vector< shared_ptr<ProtocolPartyData> >& parties){
  int nParties = _partySocket.size();
  int nActive = 0;
  vector< shared_ptr<ProtocolPartyData> > activeParties(nParties);
  for(int i=0; i<nParties; i++){
    if(_activeMask[i]){
      activeParties[nActive++] = _partySocket[i];
    }
  }
  activeParties.resize(nActive);
  parties = activeParties;
  return;
}


// thread worker function
template <class FieldType>
void BAParty<FieldType>::
exchangeBitWorker(bool sendBit, vector<bool>& recvBits,
                  int threadId, int nThreads){

  int nParties = _partySocket.size();
  char sendByte = sendBit ? 1 : 0;
  vector<char> recvBytes(nParties);

  // communicate with other parties
  for(int i=threadId; i<nParties; i+=nThreads){
    if(!_activeMask[i]){
      // skip inactive (eliminated) parties
      continue;
    }
    
    if( _myId < _partySocket[i]->getID() ){
      // write before read
      _partySocket[i]->getChannel()->write((byte*) &sendByte, 1 );
      _partySocket[i]->getChannel()->read((byte*) &(recvBytes[i]), 1 );
    }else{
      // read before write
      _partySocket[i]->getChannel()->read((byte*) &(recvBytes[i]), 1 );
      _partySocket[i]->getChannel()->write((byte*) &sendByte, 1 );
    }
  }

  // convert messages into bits
  for(int i=threadId; i<nParties; i+=nThreads){
    if(!_activeMask[i]){
      continue;
    }
    recvBits[i] = (recvBytes[i] == 1);
  }
  
  return;
}


template <class FieldType>
void BAParty<FieldType>::
exchangeBit(bool sendBit, vector<bool>& recvBits){

  vector<thread> threads(_nThread);
  for(int i=0; i<_nThread; i++){
    threads[i] = thread(&BAParty::exchangeBitWorker, this,
                        sendBit, ref(recvBits), i, _nThread);
  }

  for(int i=0; i<_nThread; i++){
    threads[i].join();
  }
  return;
}


// from BGP92: consensus() only on a bit (happiness).
// -- TODO: replace with optimized result. (``Phase King'' for now)
// -- Repeat t+1 times, each time pick a different player as king.
// -- Every player sends its bit b
// -- Every player counts the number of 1s and 0s received
// -- Every player sends two bits: (count_1 >= n-t), and (count_0 >= n-t)
// -- Every player counts the number of ``true''s for 1 and 0.
// -- Every player set its bit to (count_true_1 > t)
// -- King sends to all: its bit.
// -- Every player set its bit to King's, except when (count_true_b > n-t)
template <class FieldType>
bool BAParty<FieldType>::
consensus(bool b){
  
  int nRounds = _smallT + 1;
  int nParties = _partySocket.size();
  for(int i=0; i<nRounds; i++){
    vector<bool> receivedBits(nParties);
    // first universal excange round.
    exchangeBit(b, receivedBits);

    cout << "my (" << _myId << ") view: " << endl;
    for(int j=0; j<nParties; j++){
      cout << "p" << _partySocket[j]->getID()
           << " isActive == " << _activeMask[j] << "; ";
      cout << "bj = " << receivedBits[j] << endl;
    }
  }
  
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
