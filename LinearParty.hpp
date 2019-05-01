#ifndef LINEARPARTY_H_
#define LINEARPARTY_H_

#include <stdlib.h>
#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
#include <libscapi/include/circuits/ArithmeticCircuit.hpp>
#include <libscapi/include/infra/Measurement.hpp>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include <libscapi/include/primitives/Mersenne.hpp>
#include "ProtocolTimer.h"
#include "ECC.h"
#include "BAParty.h"
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/infra/Common.hpp>
// #include <libscapi/include/primitives/Prg.hpp>
// #include "HashEncrypt.h"
#include <emmintrin.h>
#include <thread>

using namespace std;
using namespace std::chrono;

// Plan: adapt & test existing code by George bit by bit.
// - fault localization
// - performance testing
// - (trivial) malicious party
// - optimizations

// TODO:
// 1 check that all inactive elements are 0.

// Questions:
// 1. Error checking during protocol (e.g. timeout communication?)
// 2. SCALAR and SCALAR_ADD gate types?

template <class FieldType>
// class LinearParty : public Protocol, public HonestMajority, MultiParty{
class LinearParty : public Protocol{
  
private:
  
  // ======== protocol properties ========
  int _myId;
  int _nActiveParties;
  int _bigT;
  int _smallT;
  int _nThread;
  ArithmeticCircuit _circuit;
  TemplateField<FieldType>*_field;
  FieldType _zero;
  
  // _alpha[i] is the field element of the party with id = i.
  // Note: every party has the same _alpha.
  //  ..   Same with _inputSizes
  vector<FieldType> _alpha;
  vector<FieldType> _beta;
  HIM<FieldType> _M;
  
    // subprotocol classes
  BAParty<FieldType> _baParty;
  ECC<FieldType> _eccAlpha;
  ECC<FieldType> _eccBeta;

  // manage active/inactive parties by a mask
  vector< shared_ptr<ProtocolPartyData> > _parties;
  // _activeMask[ partyId ] == 1 if partyId is active.
  vector<bool> _activeMask;
  // has size 3. each of size at most _smallT + 1
  vector< vector<int> > _partition;

  vector<int> _inputSizes;
  vector<long> _output;

  // the following are generated at the offline phase
  int _inputOffset = 0;
  vector<FieldType> _input;
  
  int _singleRandomSharesOffset = 0;
  vector<FieldType> _singleRandomShares;

  int _singleZeroSharesOffset = 0;
  vector<FieldType> _singleZeroShares;

  vector<FieldType> _wireShares;


  // ======== private functionalities ========
  
  // ---- helper funtions ----
  inline string getArg(string argKey)
  { return this->getParser().getValueByKey(this->arguments, argKey); }
  void encodeFieldElts(vector<FieldType>& input, vector<byte>& output);
  void encodeFieldElt(FieldType& input, vector<byte>& output);
  void decodeFieldElts(vector<byte>& input, vector<FieldType>& output);
  void decodeFieldElt(vector<byte>& input, FieldType& output);
  void readInputFromFile(string fileName);
  // generate a random degree d polynomial with p(0) = s
  void randomSecretePoly(FieldType& s, int d,
                         vector<FieldType>& polynomial);
  
  void expandNShareForPartition(int partitionIdx,
                                vector<FieldType>& rShareNp,
                                vector<FieldType>& rShareT);
  // expand a whole (n'-1)-share into 3 r-shares
  void expandNShare(vector<FieldType>& rShareNp, // input. a (n'-1) share
                    vector<FieldType>& rShareT1, // output, a t-share
                    vector<FieldType>& rShareT2, // output, a t-share
                    vector<FieldType>& rShareT3); // output, a t-share
  void extractNShare(vector<FieldType>& rShareT1, // input, a t-share
                     vector<FieldType>& rShareT2, // input, a t-share
                     vector<FieldType>& rShareT3, // input, a t-share
                     vector<FieldType>& rShareNp); // output, a (n'-1)-share
  void eliminateParties(const vector<int>& elimIds);
  
  // check if points lie on a degree d polynomial
  bool checkPolynomial(vector<FieldType>& y, int d, vector<FieldType>& g);
  

  // ---- helper round functions ----
  void spreadElms(vector<FieldType>& buffer,
                 int nElements, int rootId);

  void scatterForDealers(vector<FieldType>& sendElms,
                         vector<FieldType>& recvElms,
                         vector<int>& dealerIds);
  void scatterForDealersMulti(vector<vector<FieldType>>& sendElms,
                  vector<vector<FieldType>>& recvElms,
                  vector<int>& dealerIds);
  void gatherForDealers(vector<FieldType>& sendElms,
                        vector<FieldType>& recvElms,
                        vector<int>& dealerIds);
  void gatherForDealersMulti(vector<vector<FieldType>>& sendElms,
                 vector<vector<FieldType>>& recvElms,
                 vector<int>& dealerIds);
  void scatterForAll(vector<FieldType>& sendElms,
                     vector<FieldType>& recvElms);

  // ---- initialization funtions ----
  void makeField();
  void makeParties();
  void makeAlpha();
  void makeBeta();
  void makePartition();

  // ---- Batch broadcasting (BTH08 Appendix) ----
  // expand polynomial to n evaluations. 
  void prepareExpandedElms(vector<FieldType>& polynomial,
                           vector<FieldType>& expandedElms);
  // Broadcast (each dealer (k intotal) spreads T values)
  void robustBatchBroadcast(vector<FieldType>& sendElms,
                            vector< vector<FieldType> >& recvElms,
                            int nElements,  vector<int>& dealerIds);
 
  // ---- creating Random shares ----
  // create T random shares among parties. return false if fails
  // -- For a verifier party in the protocol,
  // -- ``checkVal'' is the reconstructed polinomial evaluated at 0.
  // -- For others, just set to 0.
  bool singleShareRandom(int d, vector<FieldType>& shares,
                         FieldType& checkVal);
  // instead of choosing a random to share, each party takes a input to share
  bool singleShareSecrete(FieldType& s, int d, vector<FieldType>& shares,
                          FieldType& checkVal);

  void singleShareSecreteVerify(FieldType& s, int d, 
                  vector<FieldType>& sendShares, vector<FieldType>& recvShares, vector<FieldType>& verifyShares, vector<FieldType>& toVerifyShares);

  bool singleShareSecreteVerifyOne(FieldType& s, int d, 
  vector<FieldType>& claimSendShares, vector<FieldType>& recvShares,
  vector<FieldType>& claimVerifyShares);

  void partyDispute(int party1, int party2, FieldType& sent, FieldType& received, int msgIdx,
  vector<vector<FieldType>>& mySend, vector<vector<FieldType>>& myRecv);

  // similarly create T random multiple-shares among parties
  bool multipleShareRandom(const vector<int>& degrees,
                           vector< vector<FieldType> >& shares);
  // similarly create T random shares of 0 among parties
  bool singleShareZero(int d, vector<FieldType>& shares);

  // ---- Reconstructions ----
  // reconstruct an element towards root. return 0 if fails
  FieldType reconstructPrivate(FieldType& share, int degree, int root);
  // reconstruct (at most) T elements towards all. return false if fails.
  bool reconstructPublic(vector<FieldType>& shares,
                         vector<FieldType>& results, int degree);


  // ---- creating multiplicative tuples ----
  // generate _bigT random multiplicative tuples or return false
  bool generateTuples(vector< vector<FieldType> >& a,
                      vector< vector<FieldType> >& b,
                      vector<FieldType>& c);
  
  // ---- fault localization for different subprotocols ----
  // fault localization for checkConsistency()
  void faultLocalizeConsistency(vector<int>& elimIds);
  void faultLocalizeEvalSeg(vector<int>& elimIds);

  // ---- PE versions ----
  bool generateTuplesPE(vector< vector<FieldType> >& a,
                        vector< vector<FieldType> >& b,
                        vector<FieldType> c, vector<int>& elimIds);
  

  // ---- main subprotocols ----
  
  // From new paper F_input() and F_rand()
  // -- runs linearInput() until all inputs are shared
  // -- just plug in random t-shares from offline phase input random gates
  void InputPhase();
  // From BTH08 Appendix: batch input sharing (K dealers each T inputs)
  // batch share nInputs # of input.
  // -- only works when nInputs <= _bigT
  void linearInput(vector< vector<FieldType> >& recvInputs,
                   int nInputs,  vector<int>& dealerIds);
  

  // From new paper =eval()=
  void EvalPhase();
  // computes a mult gate, consuming a multiplicative tuple
  FieldType evalMultGate(int kingId, FieldType& xShare, FieldType& yShare,
                         FieldType& aShareT, FieldType& aShareNp,
                         FieldType& bShareT, FieldType& bShareNp,
                         FieldType& cShareT, FieldType& d, FieldType& e);
  bool faultDetection(bool happiness);
  bool checkConsistency(int kingId, vector<FieldType>& vals,
                        vector<int>& elimIds);
  bool evalSeg(vector<TGate>& circuitSeg, vector<int>& elimIds);
  // generate _bigT random 4-consistent tuples or return false
  // Called QuadrupleShareRandom in the new paper
  bool generate4Tuples(vector< vector<FieldType> >& rTuple);
  bool check4Consistency();


  // From new paper =output()=
  // -- run singleShareZero() for each batch of (_bigT) output gates.
  // -- and then reconstruct output share to the owner
  void OutputPhase();
  
public:

  LinearParty(int argc, char* argv[]);
  ~LinearParty();

  // ======== inherited functionalites ========
  bool hasOffline() override {return true;}
  bool hasOnline() override {return true;}
  void run() override; // <- runs =main()= from the new paper
  void runOffline() override;
  void runOnline() override;

  // ======== testing functionalities ========

  
    
};


static inline void 
fillInputSizes(int myId, int nParties, ArithmeticCircuit& circuit, // input
               vector<int>& sizes){ // output
  int nInputGates = circuit.getNrOfInputGates();
  const vector<TGate> gates = circuit.getGates();

  sizes.clear();
  sizes.resize(nParties);
  int i=0;
  while(i < nInputGates){
    if(gates[i].gateType == INPUT){
      sizes[gates[i++].party]++;
    }
  }
  
  return;
}

static inline void 
lastNActiveParties(int lastN, const vector<bool>& mask,
                   vector<int>& activeIds) {
  int nPartiesInc = mask.size();
  int lastId = nPartiesInc-1;
  activeIds.resize(lastN);
  while (lastN > 0 && lastId >= 0) {
    if (mask[lastId]) {
      activeIds[lastN-1] = lastId;
      lastN--;
    }
    lastId--;
  }
  return;
}

static inline void
firstNActiveParties(int firstN, const vector<bool>& mask,
                    vector<int>& activeIds) {
  int firstId = 0;
  int idx = 0;
  activeIds.resize(firstN);
  while (idx < firstN) {
    if (mask[firstId]) {
      activeIds[idx] = firstId;
      idx++;
    }
    firstId++;
  }
  return;
}

template <class FieldType>
static inline void
firstNActiveValues(int firstN, const vector<bool> &mask,
                   vector<FieldType> &totalValues,
                   vector<FieldType> &activeValues) {
  
  vector<int> activeIds;
  firstNActiveParties(firstN, mask, activeIds);
  activeValues.resize(firstN);
  for (int i = 0; i < firstN; i++) {
    activeValues[i] = totalValues[ activeIds[i] ];
  }
  return;
}

template <class FieldType>
static inline void
lastNActiveValues(int lastN, const vector<bool> &mask,
                  vector<FieldType> &totalValues,
                  vector<FieldType> &activeValues) {
  
  vector<int> activeIds;
  lastNActiveParties(lastN, mask, activeIds);
  activeValues.resize(lastN);
  for (int i = 0; i < lastN; i++) {
    activeValues[i] = totalValues[ activeIds[i] ];
  }
  return;
}

template<class FieldType>
void LinearParty<FieldType>::
makeField(){
  string fieldType = getArg("fieldType");
  if(fieldType == "ZpMersenne31"){
    _field = new TemplateField<FieldType>(2147483647);
  }else if(fieldType == "ZpMersenne61"){
    _field = new TemplateField<FieldType>(0);
  }else if(fieldType == "ZpKaratsuba"){
    _field = new TemplateField<FieldType>(0);
  }else if(fieldType == "GF2m"){
    _field = new TemplateField<FieldType>(8);
  }else if(fieldType == "Zp"){
    _field = new TemplateField<FieldType>(2147483647);
  }
  
  _zero = *(_field->GetZero());
  return;  
}


template<class FieldType>
void LinearParty<FieldType>::
makeParties(){

  MPCCommunication comm;
  boost::asio::io_service io_service;
  _parties =
    comm.setCommunication(io_service, _myId, _nActiveParties,
                          getArg("partiesFile"));

  _activeMask.clear();
  _activeMask.resize(_nActiveParties, true);
  return;
}

template<class FieldType>
void LinearParty<FieldType>::
makeAlpha(){
  _alpha.resize(_nActiveParties);
  for(int i = 0; i < _nActiveParties; i++) {
    _alpha[i] = _field->GetElement(i + 1);
  }
  return;
}


template<class FieldType>
void LinearParty<FieldType>::
makeBeta(){
  _beta.resize(_nActiveParties);
  for(int i = 0; i < _nActiveParties; i++) {
    _beta[i] = _field->GetElement(i + 1 + _nActiveParties);
  }
  return;
}

template<class FieldType>
void LinearParty<FieldType>::
makePartition(){
  // divide all party ids into 3 partitions, each of size
  // at most _smallT + 1;
  _partition.resize(3);
  int curId = 0;
  int maxSize = _smallT + 1;
  int partiesLeft =  _parties.size() + 1;
  int curPartitionSize = 0;

  // just assign partitions in id order.
  for (int j = 0; j < 3; j++) {
    partiesLeft -= curPartitionSize;
    curPartitionSize = maxSize <= partiesLeft ? maxSize : partiesLeft;
    _partition[j].resize(curPartitionSize);
    for (int i = 0; i < curPartitionSize; i++) {
      _partition[j][i] = curId++;
    }
  }

  return;
}

template <class FieldType>
LinearParty<FieldType>::
LinearParty(int argc, char* argv[])
  : Protocol("MPCLinearComm", argc, argv)
{
  // fill basic class properties
  _myId = stoi(getArg("partyID"));
  _nActiveParties = stoi(getArg("partiesNumber"));
  _nThread = stoi(getArg("numThreads"));
  _smallT = (_nActiveParties -1)/ 3;
  _bigT = _nActiveParties - 2*_smallT;
  
  makeField();
  makeParties();
  makeAlpha();
  makeBeta();
  makePartition();

  // initialize subprotocol classes
  _M.allocate(_nActiveParties, _nActiveParties, _field);
  _M.InitHIMByVectors(_alpha, _beta);
  _baParty.setParties(_parties, _myId);
  _baParty.setNumThreads(_nThread);
  _baParty.setAlphaBeta(_alpha, _beta);
  _baParty.setHIM(_M);
  _eccAlpha.setField(_field); // need to set field first!
  _eccAlpha.setAlpha(_alpha);
  _eccBeta.setField(_field);
  _eccBeta.setAlpha(_beta);
  
      
  // build _circuit
  _circuit.readCircuit( (getArg("circuitFile")).c_str() );
  _circuit.reArrangeCircuit();
  _input.resize(_circuit.getNrOfInputGates());
  _output.resize(_circuit.getNrOfOutputGates());
  
  // build _inputSizes
  fillInputSizes(_myId, _parties.size() + 1, _circuit, _inputSizes);
  readInputFromFile(getArg("inputFile"));
  return;
}

template <class FieldType>
LinearParty<FieldType>::
~LinearParty(){
  delete _field;
}


template <class FieldType>
void LinearParty<FieldType>::
run(){

  auto start = high_resolution_clock::now();
  
  runOffline();
  runOnline();

  auto end = high_resolution_clock::now();
  auto duration =
    duration_cast<milliseconds>(end-start).count();

  cout << "time in milliseconds initializationPhase: "
       << duration << endl;

  return;
}


template <class FieldType>
void LinearParty<FieldType>::
encodeFieldElts(vector<FieldType>& input, vector<byte>& output){
  int fieldByteSize = _field->getElementSizeInBytes();
  output.resize(fieldByteSize * input.size());
  for(int j = 0; j < input.size(); j++) {
    _field->elementToBytes(output.data() + (j * fieldByteSize), input[j]);
  }
}

template <class FieldType>
void LinearParty<FieldType>::
encodeFieldElt(FieldType& input, vector<byte>& output){
  int fieldByteSize = _field->getElementSizeInBytes();
  output.resize(fieldByteSize);
  _field->elementToBytes(output.data(), input);
}

template <class FieldType>
void LinearParty<FieldType>::
decodeFieldElts(vector<byte>& input, vector<FieldType>& output){
  int fieldByteSize = _field->getElementSizeInBytes();
  // if(input.size() % fieldByteSize) {
  //   cout << "input byte misalignment" << endl;
  // }
  output.resize(input.size() / fieldByteSize);

  for(int j = 0; j < output.size(); j++) {
    output[j] = _field->bytesToElement(input.data() + fieldByteSize * j);
  }
}

template <class FieldType>
void LinearParty<FieldType>::
decodeFieldElt(vector<byte>& input, FieldType& output){
  int fieldByteSize = _field->getElementSizeInBytes();
  output = _field->bytesToElement(input.data());
}

template<class FieldType>
void LinearParty<FieldType>::
readInputFromFile(string fileName){

  ifstream inputFile(fileName);
  int inputSize = _input.size();
  long curInput;

  // read available inputs to fill _input
  for(int i=0; i< inputSize && !(inputFile.eof()); i++){
    inputFile >> curInput;
    _input[i] = _field->GetElement(curInput);
  }
  return;
}


// spread from a root to all
// a specified number of field elements
template <class FieldType>
void LinearParty<FieldType>::
spreadElms(vector<FieldType>& buffer,
           int nElements, int rootId){

  vector<byte> msg;
  int nParties = _parties.size();
  int msgSize = _field->getElementSizeInBytes() * nElements;

  if (rootId == _myId) {
    // I'm the sender
    encodeFieldElts(buffer, msg);
    msg.resize(msgSize);
    _baParty.spreadMsg(msg, msgSize, rootId);
  } else {
    // I'm a receiver
    _baParty.spreadMsg(msg, msgSize, rootId);
    decodeFieldElts(msg, buffer);
  }
  return;
}


template <class FieldType>
void LinearParty<FieldType>::
prepareExpandedElms(vector<FieldType>& polynomial,
                    vector<FieldType>& expandedElms){

  int nPartiesInc = _parties.size()+1;
  expandedElms.resize(nPartiesInc);
  for (int i = 0; i < nPartiesInc; i++) {
    expandedElms[i] = _eccAlpha.evalPolynomial(_alpha[i], polynomial);
  }

  return;
}


// Broadcast (each dealer (k intotal) spreads T values)
// -- Every dealer expand T values into n by interpolating
// -- Every dealer distribute 1 value to each party
// -- Run BroadcastForP() for k values (1 from each dealer)
// -- Every player reconstruct T values for each dealer w/ ECC
// -- note that we must have T > n - 2t.
template <class FieldType>
void LinearParty<FieldType>::
robustBatchBroadcast(vector<FieldType>& sendElms, // input
                     vector< vector<FieldType> >& recvElms, // output
                     int nElements,  vector<int>& dealerIds){

  // sort dealerIds. 
  int nParties = _parties.size();
  int nDealers = dealerIds.size();
  int elmSize = _field->getElementSizeInBytes();

  // get one element from each dealer Id
  vector<FieldType> sendExpandedElms;
  prepareExpandedElms( sendElms, sendExpandedElms );
  vector<FieldType> recvExpandedElms(nDealers);
  scatterForDealers(sendExpandedElms, recvExpandedElms, dealerIds);
  
  // each party broadcast all received shares
  int msgSize = elmSize * nDealers;
  vector<byte> sendMsg;
  vector< vector<byte> > recvMsgs;
  encodeFieldElts(recvExpandedElms, sendMsg);
  _baParty.broadcastMsgForAll(sendMsg, recvMsgs, msgSize);
  
  // each party decode received shares
  vector< vector<FieldType> > // each row is an expanded code.
    allExpandedElms(nDealers, vector<FieldType>(nParties+1));
  for (int i = 0; i < nDealers; i++) {
    allExpandedElms[i][ _myId ] = recvExpandedElms[i];
  }

  vector<FieldType> tmpVector;
  for (int i = 0; i < nParties; i++) {
    int partyId = _parties[i]->getID();
    decodeFieldElts(recvMsgs[i], tmpVector);
    for (int j = 0; j < nDealers; j++) {
      allExpandedElms[j][partyId] = tmpVector[j];
    }
  }

  // use ECC to reconstruct each dealer's elms
  recvElms.resize(nDealers);
  for (int i = 0; i < nDealers; i++) {
    _eccAlpha.reconstruct(allExpandedElms[i], nElements-1, recvElms[i]);
  }

  return;
}

template <class FieldType>
void LinearParty<FieldType>::
randomSecretePoly(FieldType &s, int degree,
                  vector<FieldType> &polynomial) {
  polynomial.resize(degree+1);
  polynomial[0] = s;
  for (int i = 0; i < degree; i++) {
    polynomial[i+1] = _field->Random();
  }

  return;
}

template <class FieldType>
void LinearParty<FieldType>::
scatterForDealers(vector<FieldType>& sendElms,
                  vector<FieldType>& recvElms,
                  vector<int>& dealerIds) {
  int nDealers = dealerIds.size();
  int nParties = _parties.size();
  int nPartiesInc = _parties.size() +1;
  sort(dealerIds.begin(), dealerIds.end());
  recvElms.resize(nDealers);
  sendElms.resize(nPartiesInc);
   
  int msgSize = _field->getElementSizeInBytes();
  vector<byte> recvMsg(msgSize);
  vector< vector<byte> > sendMsgs(nParties, vector<byte>(msgSize));

  for(int i=0; i<nDealers; i++){
    int dealerId = dealerIds[i];

    if ( dealerId == _myId) {
      recvElms[i] = sendElms[_myId];
      for (int j = 0; j < nParties; j++) {
        encodeFieldElt(sendElms[ _parties[j]->getID() ], sendMsgs[j]);
      }
      _baParty.scatterMsg(sendMsgs, recvMsg, msgSize, dealerId);
      
    } else {
      // receive message and store
      _baParty.scatterMsg(sendMsgs, recvMsg, msgSize, dealerId);
      decodeFieldElt(recvMsg, recvElms[i]);
    }
  }

  return;
}

template <class FieldType>
void LinearParty<FieldType>::
scatterForDealersMulti(vector<vector<FieldType>>& sendElms,
                  vector<vector<FieldType>>& recvElms,
                  vector<int>& dealerIds) {
  int nDealers = dealerIds.size();
  int nParties = _parties.size();
  int nPartiesInc = _parties.size() +1;
  sort(dealerIds.begin(), dealerIds.end());
  recvElms.resize(nDealers);
  sendElms.resize(nPartiesInc);
   
  int msgSize = _field->getElementSizeInBytes() * sendElms[0].size();
  vector<byte> recvMsg(msgSize);
  vector< vector<byte> > sendMsgs(nParties, vector<byte>(msgSize));

  for(int i=0; i<nDealers; i++){
    int dealerId = dealerIds[i];

    if ( dealerId == _myId) {
      recvElms[i] = sendElms[_myId];
      for (int j = 0; j < nParties; j++) {
        encodeFieldElts(sendElms[ _parties[j]->getID() ], sendMsgs[j]);
      }
      _baParty.scatterMsg(sendMsgs, recvMsg, msgSize, dealerId);
      
    } else {
      // receive message and store
      _baParty.scatterMsg(sendMsgs, recvMsg, msgSize, dealerId);
      decodeFieldElts(recvMsg, recvElms[i]);
    }
  }

  return;
}

template <class FieldType>
void LinearParty<FieldType>::
gatherForDealers(vector<FieldType>& sendElms,
                 vector<FieldType>& recvElms,
                 vector<int>& dealerIds) {
  int nDealers = dealerIds.size();
  int nParties = _parties.size();
  int nPartiesInc = _parties.size() +1;

  
  sendElms.resize(nDealers);
  recvElms.clear(); // <---- non-dealer can check recvElms.size()
  // <---------------------- to know that they are not dealers
  
  int msgSize = _field->getElementSizeInBytes();
  vector<byte> sendMsg(msgSize);
  vector< vector<byte> > recvMsgs;
  for(int i=0; i<nDealers; i++){
    int dealerId = dealerIds[i];

    if ( dealerId == _myId) {
      recvElms.resize(nPartiesInc);
      recvElms[_myId] = sendElms[i];
      _baParty.gatherMsg(sendMsg, recvMsgs, msgSize, dealerId);
      
      for (int j = 0; j < nParties; j++) {
        decodeFieldElt(recvMsgs[j], recvElms[ _parties[j]->getID() ]);
      }
      
    } else {
      encodeFieldElt(sendElms[i], sendMsg);
      _baParty.gatherMsg(sendMsg, recvMsgs, msgSize, dealerId);
    }
  }

  return;
}

template <class FieldType>
void LinearParty<FieldType>::
gatherForDealersMulti(vector<vector<FieldType>>& sendElms,
                 vector<vector<FieldType>>& recvElms,
                 vector<int>& dealerIds) {
  int nDealers = dealerIds.size();
  int nParties = _parties.size();
  int nPartiesInc = _parties.size() +1;

  
  sendElms.resize(nDealers);
  recvElms.clear(); // <---- non-dealer can check recvElms.size()
  // <---------------------- to know that they are not dealers
  
  int msgSize = _field->getElementSizeInBytes() * sendElms[0].size();
  vector<byte> sendMsg(msgSize);
  vector< vector<byte> > recvMsgs;
  for(int i=0; i<nDealers; i++){
    int dealerId = dealerIds[i];

    if ( dealerId == _myId) {
      recvElms.resize(nPartiesInc);
      recvElms[_myId] = sendElms[i];
      _baParty.gatherMsg(sendMsg, recvMsgs, msgSize, dealerId);
      
      for (int j = 0; j < nParties; j++) {
        decodeFieldElts(recvMsgs[j], recvElms[ _parties[j]->getID() ]);
      }
      
    } else {
      encodeFieldElts(sendElms[i], sendMsg);
      _baParty.gatherMsg(sendMsg, recvMsgs, msgSize, dealerId);
    }
  }

  return;
}


template <class FieldType>
void LinearParty<FieldType>::
scatterForAll(vector<FieldType>& sendElms,
              vector<FieldType>& recvElms) {

  int nPartiesInc = _parties.size() + 1;
  vector<int> dealers;
  for (int i = 0; i < nPartiesInc; i++) {
    if (_activeMask[i]) {
      dealers.push_back(i);
    }
  }

  scatterForDealers(sendElms, recvElms, dealers);
  return;
}

template <class FieldType>
bool LinearParty<FieldType>::
checkPolynomial(vector<FieldType>& y, int degree, vector<FieldType>& g) {

  bool happiness = true;
  int nPoints = y.size();
  // assert(y.size() == _parties().size()+1)
  g.clear();
  vector<FieldType> xTmp(degree+1);
  vector<FieldType> yTmp(degree+1);

  firstNActiveValues(degree+1, _activeMask, _alpha, xTmp);
  firstNActiveValues(degree+1, _activeMask, y, yTmp);
  _eccAlpha.interpolate(xTmp, yTmp, g);
  // check the rest
  for (int i = 0; i < nPoints; i++) {
    if (_activeMask[i]) {
      happiness &=
        (_eccAlpha.evalPolynomial(_alpha[i], g) == y[i]);
    }
  }
  
  return happiness;
}


// create T random shares among parties. return false if fails
// Note: non-robust.
// -- every party takes input a secrete, and share with all others.
// -- every party applies HIM on received shares
// -- every party send transformedShares of T+1 ... n' to crrspnd. parties
// -- the (n'-T) verifiers verify (n'-T) received shares
//TODO: Player elimination version
template <class FieldType>
bool LinearParty<FieldType>::
singleShareSecrete(FieldType& s, int d, vector<FieldType>& shares,
                  FieldType& checkVal) {

  bool happiness = true;
  int nPartiesInc = _parties.size() + 1;

  // -- every party choose a random secrete
  // FieldType r = _field->Random();
  FieldType r = s;
  vector<FieldType> f;
  vector<FieldType> sendShares(nPartiesInc);


  randomSecretePoly(r, d, f); //Creates Random Bits in F


  for (int i = 0; i < nPartiesInc; i++) {
    sendShares[i] = _eccAlpha.evalPolynomial(_alpha[i], f);
  }
  // -- and share with all others.
  vector<FieldType> recvShares(nPartiesInc);


  scatterForAll(sendShares, recvShares); //communication


  // -- every party applies HIM on received shares
  vector<FieldType> transformedShares(nPartiesInc);
  _M.MatrixMult(recvShares, transformedShares);

  // -- every party send transformedShares of T+1 ... n' to
  // -- the corresponding parties
  int nVerifiers = _nActiveParties - _bigT;
  vector<int> dealers(nVerifiers);
  lastNActiveParties(nVerifiers, _activeMask, dealers);
  vector<FieldType> verifyShares(nVerifiers);
  lastNActiveValues(nVerifiers, _activeMask, transformedShares, verifyShares);
  vector<FieldType> toVerifyShares; // received if I'm a verifier


  gatherForDealers(verifyShares, toVerifyShares, dealers); //communication


  // -- the (n'-T) verifiers verify (n'-T) received shares
  checkVal = *(_field->GetZero()); // default value
  if (toVerifyShares.size() > 0) {
    vector<FieldType> g;
    happiness &= checkPolynomial(toVerifyShares, d, g);
    
    // compute the checkVal
    if (happiness) {
      checkVal = _eccAlpha.evalPolynomial(*(_field->GetZero()), g);
    }
  }

  if(!faultDetection(happiness)) { //byzantine unhappy, verify consistency and sharing
    singleShareSecreteVerify(s, d, sendShares, recvShares, verifyShares, toVerifyShares);
  }

  // -- output the first _bigT shares from active parties
  firstNActiveValues(_bigT, _activeMask, transformedShares, shares);
  return happiness;
}

template <class FieldType> 
void LinearParty<FieldType>::
singleShareSecreteVerify(FieldType& s, int d, 
                  vector<FieldType>& sendShares, vector<FieldType>& recvShares, 
                  vector<FieldType>& verifyShares, vector<FieldType>& toVerifyShares) {

  FieldType x; //just some filler args
  vector<int> judge;
  firstNActiveParties(1, _activeMask, judge);

  int nVerifiers = _nActiveParties - _bigT;
  vector<int> dealers(nVerifiers);
  lastNActiveParties(nVerifiers, _activeMask, dealers);

  vector<FieldType> random1;
  random1.push_back(s);
  vector<FieldType> randomElts;

  vector<vector<FieldType>> sendSharesAll;
  vector<vector<FieldType>> recvSharesAll;
  vector<vector<FieldType>> verifySharesAll;
  vector<vector<FieldType>> toVerifySharesAll;

  vector<vector<FieldType>> sendShares1;
  sendShares1.push_back(sendShares);
  
  vector<vector<FieldType>> recvShares1;
  recvShares1.push_back(recvShares);
  
  vector<vector<FieldType>> verifyShares1;
  vector<FieldType> verifySharesComplete(_parties.size() + 1, *(_field->GetZero()));
  for(unsigned i = 0; i < dealers.size(); i++) {
    verifySharesComplete[dealers[i]] = verifyShares[i];
  }
  verifyShares1.push_back(verifyShares);
  
  vector<vector<FieldType>> toVerifyShares1;
  toVerifyShares.resize(_parties.size() + 1); //so consistent size with all dealers
  toVerifyShares1.push_back(toVerifyShares);

  //send to the judge - sendShares, recvShares, verifyShares, toVerifyShares(if dealer)
  gatherForDealers(random1, randomElts, judge);
  gatherForDealersMulti(sendShares1, sendSharesAll, judge); //everyone sending their data to judge. 
  gatherForDealersMulti(recvShares1, recvSharesAll, judge); 
  gatherForDealersMulti(verifyShares1, verifySharesAll, judge); 
  gatherForDealersMulti(toVerifyShares1, toVerifySharesAll, judge); 
  //toVerifyShares may be different sizes for dealers / non dealers

  vector<vector<FieldType>> sent;
  sent.push_back(sendShares);
  sent.push_back(verifyShares);

  vector<vector<FieldType>> recv;
  recv.push_back(recvShares);
  recv.push_back(toVerifyShares);

  if(_myId == judge[0]) { //I am the judge
    //checkConsistency
    for(int i = 0; i < sendShares.size(); i++) {
      for(int j = 0; j < recvShares.size(); j++) {
        if(sendSharesAll[i][j] != recvSharesAll[j][i]) {
          //Dispute!!
          partyDispute(i,j,sendSharesAll[i][j], recvSharesAll[j][i], 0, sent, recv);
          return;
        }
      }
    }

    for(int i = 0; i < verifyShares.size(); i++) {
      for(int j = 0; j < toVerifyShares.size(); j++) {
        if(verifySharesAll[i][j] != toVerifySharesAll[j][i]) {
          partyDispute(i,j,verifySharesAll[i][j], toVerifySharesAll[j][i], 1, sent, recv);
          return;
        }
      }
    }
    
    //simulate computation
    for(int i = 0; i < _parties.size() + 1; i++) {
      if(_activeMask[i]) {
        if(!singleShareSecreteVerifyOne(randomElts[i], d, sendSharesAll[i], recvSharesAll[i], verifySharesAll[i])) {
          partyDispute(_myId, i, x, x, 0, sent, recv);
        }
      }
    }


  } else { //listen for disputes
    partyDispute(0,0,x,x,0, sent, recv);
  }
}


template <class FieldType>
bool LinearParty<FieldType>::
singleShareSecreteVerifyOne(FieldType& s, int d, 
  vector<FieldType>& claimSendShares, vector<FieldType>& recvShares,
  vector<FieldType>& claimVerifyShares) {
  //first check claimed sendShares are on degree d polynomial
  vector<FieldType> poly;
  _eccAlpha.interpolate(claimSendShares, poly);
  if(poly.size() > d + 1) {
    return false;
  }
  //confirm rec generates claimed verify shares.

  vector<FieldType> transformedShares(_parties.size() + 1);
  _M.MatrixMult(recvShares, transformedShares);

  // -- every party send transformedShares of T+1 ... n' to
  // -- the corresponding parties
  int nVerifiers = _nActiveParties - _bigT;
  vector<int> dealers(nVerifiers);
  lastNActiveParties(nVerifiers, _activeMask, dealers);
  vector<FieldType> verifyShares(nVerifiers);
  lastNActiveValues(nVerifiers, _activeMask, transformedShares, verifyShares);
  
  if(verifyShares.size() != claimVerifyShares.size()) {
    return false;
  }
  for(int i = 0; i < verifyShares.size(); i++) {
    if(verifyShares[i] != claimVerifyShares[i]) {
      return false;
    }
  }

  return true;
}

template <class FieldType>
void LinearParty<FieldType>::
partyDispute(int party1, int party2, FieldType& sent, FieldType& received, int msgIdx,
  vector<vector<FieldType>>& mySend, vector<vector<FieldType>>& myRecv) {
  vector<int> judge;
  firstNActiveParties(1, _activeMask, judge);

  //assumes < 255 parties
  byte party1b = (byte)party1;
  byte party2b = (byte)party2;
  byte msgIdxb = (byte)msgIdx;
  vector<byte> msg = { party1b, party2b, msgIdxb };
  vector<vector<byte>> msgs(_parties.size(), msg);
  vector<byte> recMsg;

  _baParty.scatterMsg(msgs, recMsg, 2, judge[0]);
  vector<FieldType> log;
  log.push_back(sent);
  log.push_back(received);
  vector<vector<FieldType>> msg2(_parties.size() + 1, log);
  vector<vector<FieldType>> recMsg2;
  scatterForDealersMulti(msg2, recMsg2, judge);

  party1 = (int)msg[0];
  party2 = (int)msg[1];
  msgIdx = (int)msg[2];

  sent = recMsg2[0][0];
  received = recMsg2[0][1];

  if(party1 == judge[0] || party2 == judge[0]) {
    vector<int> eliminated;
    eliminated.push_back(party1);
    eliminated.push_back(party2);
    eliminateParties(eliminated);
  } else {
    FieldType dispute = *(_field->GetZero()); //lazy, this is just a boolean
    vector<FieldType> disputeAll;
    if(_myId == party1) {
      if(mySend[msgIdx][party2] != sent) {
        dispute = *(_field->GetOne());
      }
    } else if(_myId == party2) {
      if(myRecv[msgIdx][party1] != received) {
        dispute = *(_field->GetOne());
      }
    }
    vector<int> conflicters = { party1, party2 };
    vector<FieldType> disputev(2, dispute);
    scatterForDealers(disputev, disputeAll, conflicters);

    vector<int> eliminated;
    if(disputeAll[party1] != *(_field->GetZero())) {
      eliminated.push_back(judge[0]);
      eliminated.push_back(party1);
    } else if (disputeAll[party2] != *(_field->GetZero())) {
      eliminated.push_back(judge[0]);
      eliminated.push_back(party2);
    } else {
      eliminated.push_back(party1);
      eliminated.push_back(party2);
    }

    eliminateParties(eliminated);
  }
}

// choose a random secrete, and run singleShareSecrete
template <class FieldType>
bool LinearParty<FieldType>::
singleShareRandom(int d, vector<FieldType>& shares, FieldType& checkVal) {
  bool happiness = true;
  FieldType r;
  
  do {
    r = _field->Random();
  } while(!singleShareSecrete(r, d, shares, checkVal)); //if false, then parties eliminated. Rerun

  return true;
}

// similarly create T random multiple-shares among parties
// -- For a verifyer, check degree of received shares
// -- and also that the reconstruced values are equal
template <class FieldType>
bool LinearParty<FieldType>::
multipleShareRandom(const vector<int>& degrees,
                    vector< vector<FieldType> > &shares) {
  bool happiness = true;
  int nShares = degrees.size();
  shares.resize(nShares);
  FieldType r = _field->Random();
  vector<FieldType> checkVals(nShares);

  for (int i = 0; i < nShares; i++) {
    int d = degrees[i];
    happiness &= singleShareSecrete(r, d, shares[i], checkVals[i]);
  }

  for (int i = 0; i < nShares-1; i++) {
    happiness &= (checkVals[i] == checkVals[i+1]);
  }

  return happiness;
}

// similarly create T random shares of 0 among parties
template <class FieldType>
bool LinearParty<FieldType>::
singleShareZero(int d, vector<FieldType>& shares) {
  bool happiness = true;
  vector<int> degrees(2, d);
  vector< vector<FieldType> > doubleShares;
  happiness &= multipleShareRandom(degrees, doubleShares);

  shares.resize(_bigT);
  for (int i = 0; i < _bigT; i++) {
    shares[i] = doubleShares[0][i] - doubleShares[1][i];
  }
  
  return happiness;
}

// reconstruct an element towards root. return 0 if fails
// Note: robust if degree <= T.
template <class FieldType>
FieldType LinearParty<FieldType>::
reconstructPrivate(FieldType& share, int degree, int root){
  FieldType result = _zero;

  vector<byte> msg;
  vector< vector<byte> > recvMsgs;

  int msgSize = _field->getElementSizeInBytes();
  int nPartiesInc = _parties.size() + 1;
  int nParties  = _parties.size();
  if (_myId == root) {
    // root: receive shares, and decode
    vector<FieldType> shares(nPartiesInc, _zero);
    shares[_myId] = share;
    
    recvMsgs.resize(nParties);
    _baParty.gatherMsg(msg, recvMsgs, msgSize, root);

    for (int i = 0; i < nParties; i++) {
      if (_activeMask[ _parties[i]->getID() ]) {
        decodeFieldElt(recvMsgs[i], shares[ _parties[i]->getID() ]);
      }
    }

    vector<FieldType> polynomial;
    bool success = _eccAlpha.reconstruct(shares, degree, polynomial);

    // debug TODO: remove
    if (!success) {
      cout << "degree is " <<  degree << endl;
      for (auto e : shares) {
        cout << e << " ";
      }
      cout << endl;
    }

    result = _eccAlpha.evalPolynomial(_zero, polynomial);

  } else {
    // non-root party: send my share
    encodeFieldElt(share, msg);
    _baParty.gatherMsg(msg, recvMsgs, msgSize, root);
    
  }

  return result;
}
// reconstruct (at most) _bigT elements towards all. return false if fails.
// Note: robust if degree <= _smallT.
// -- treat shares as a polynomial
// -- expand to n points by evaluating at beta
// -- send to corresponding parties (all-to-all)
// -- reconstrcut f according to received, and send f(0) (all-to-all)
// -- reconstruct g according to received,
// -- and coefficients are reconstruced values.
template <class FieldType>
bool LinearParty<FieldType>::
reconstructPublic(vector<FieldType>& shares,
                  vector<FieldType>& results, int degree){

  int nParties = _parties.size();
  int nPartiesInc = _parties.size() + 1;
  int msgSize = _field->getElementSizeInBytes();
  int batchSize = shares.size();

  if (batchSize == 0) {
    return true;
  }

  // -- treat shares as a polynomial
  // -- expand to n points by evaluating at beta
  vector<FieldType> uSend(nPartiesInc);
  for (int i = 0; i < nPartiesInc; i++) {
    uSend[i] = _eccAlpha.evalPolynomial(_beta[i], shares);
  }

  // -- send to corresponding parties (all-to-all)
  vector<FieldType> uRecv(nPartiesInc);
  scatterForAll(uSend, uRecv);


  // -- reconstrcut f according to received, and send f(0) (all-to-all)
  vector<FieldType> f;
  bool success = _eccAlpha.reconstruct(uRecv, degree, f);
  if (!success) {
    return false;
  }

  FieldType f0;
  f0 = _eccAlpha.evalPolynomial(_zero, f);
  vector<byte> sendMsg(msgSize);
  vector< vector<byte> > recvMsgs(nParties, vector<byte>(msgSize));
  encodeFieldElt(f0, sendMsg);

  _baParty.broadcastMsgForAll(sendMsg, recvMsgs, msgSize);

  // -- reconstruct g according to received
  // -- and coefficients are reconstruced values.
  vector<FieldType> recvF0(nPartiesInc);
  recvF0[_myId] = f0;
  for (int i = 0; i < nParties; i++) {
    if (!_activeMask[ _parties[i]->getID() ]) {
      continue;
    }

    decodeFieldElt(recvMsgs[i], recvF0[ _parties[i]->getID() ]);
  }
  success = _eccBeta.reconstruct(recvF0, batchSize-1, results);
  return success;
}


template <class FieldType>
bool LinearParty<FieldType>::
generateTuples(vector< vector<FieldType> >& a,
               vector< vector<FieldType> >& b,
               vector<FieldType>& c){
  bool happiness = true;

  // -- create triple random shares for a, b, and r w/ appropriate degrees 
  vector< vector<FieldType> > rShares;
  int smallTp = (_nActiveParties - _bigT) / 2;
  vector<int> abDegrees{_smallT, _nActiveParties-1, smallTp};
  vector<int> rDegrees{_smallT, 2*smallTp}; // no need for (n'-1)-share
  happiness &= multipleShareRandom(abDegrees, a);
  happiness &= multipleShareRandom(abDegrees, b);
  happiness &= multipleShareRandom(rDegrees, rShares);

  // -- locally compute 2tp-shares of c := ab
  // -- locally compute 2tp-shares of d := c - r
  vector<FieldType> dShares2Tp(_bigT);
  for (int i = 0; i < _bigT; i++) {
    dShares2Tp[i] = a[2][i] * b[2][i] - rShares[1][i];
  }

  // -- reconstruct publicly ds
  vector<FieldType> reconsDs(_bigT);
  happiness &= reconstructPublic(dShares2Tp, reconsDs, 2*smallTp);

  // -- compute c's and output a, b, c
  a.resize(2); // just take the t-shares and (n'-1)-shares
  b.resize(2); // for a and b. 
  c.resize(_bigT);
  for (int i = 0; i < _bigT; i++) {
    c[i] = rShares[0][i] + reconsDs[i];
  }

  return happiness;
}


// TODO implement
template <class FieldType>
bool LinearParty<FieldType>::
generateTuplesPE(vector< vector<FieldType> >& a,
                 vector< vector<FieldType> >& b,
                 vector<FieldType> c, vector<int>& elimIds){

  // -- invoke TripleShareRandomPE. three times as before
  // ---- if fails, just return elimIds.

  // -- invode ReconsPE. as before
  // ---- if fails, just return elimIds.

  // -- also remember to record all info needed for check-4-consistency
  cout << "generateTuplesPE not implemented yet"  << endl;
  assert(false);
  return true;
}


// computes a mult gate, consuming a multiplicative tuple
// -- all party computes (n-1)-share d_np := x_t + a_np
// -- all party computes (n-1)-share e_np := y_t + b_np
// -- all parties send shares of d and e to king
// -- king sends back reconstructed d and e
// -- record d and e distributed by kingId
// -- all parties compute z_t := d*e - d*b_t - e*a_t + c_t
template <class FieldType>
FieldType LinearParty<FieldType>::
evalMultGate(int kingId, FieldType &xShare, FieldType &yShare,
             FieldType &aShareT, FieldType &aShareNp,
             FieldType &bShareT, FieldType &bShareNp,
             FieldType &cShareT, FieldType &d, FieldType &e ) {

  // -- all party computes (n-1)-share d_np := x_t + a_np
  // -- all party computes (n-1)-share e_np := y_t + b_np
  FieldType dShareNp = xShare + aShareNp;
  FieldType eShareNp = yShare + bShareNp;

  // -- all parties send shares of d and e to king
  // -- king sends back reconstructed d and e
  vector<int> dealer{kingId};
  vector<FieldType> sendElm(1);
  vector<FieldType> recvElms;
  vector<FieldType> reconsDE(2); // [0] is d, and [1] is e
  vector<FieldType> g;

  sendElm[0] = dShareNp;
  gatherForDealers(sendElm, recvElms, dealer);
  if (kingId == _myId) {
    _eccAlpha.interpolate(recvElms, g);
    reconsDE[0] = _eccAlpha.evalPolynomial(*(_field->GetZero()), g);
  }
  sendElm[0] = eShareNp;
  gatherForDealers(sendElm, recvElms, dealer);
  if (kingId == _myId) {
    _eccAlpha.interpolate(recvElms, g);
    reconsDE[1] = _eccAlpha.evalPolynomial(*(_field->GetZero()), g);
  }
  
  spreadElms(reconsDE, 2, kingId);
  // -- record d and e distributed by kingId
  d = reconsDE[0];
  e = reconsDE[1];

  // -- all parties compute z_t := d*e - d*b_t - e*a_t + c_t
  FieldType zShare = reconsDE[0] * reconsDE[1] -
    reconsDE[0] * bShareT - reconsDE[1] * aShareT + cShareT;
  return zShare;
}

// detect that any party is unhappy
// -- each party sends happiness to all others
// -- set to unhappy if sees any unhappy
// -- all parties run consensus on happiness
// -- NOTE: assuming broadcast channel. no need for consensus
template <class FieldType>
bool LinearParty<FieldType>::
faultDetection(bool happiness) {
  int nParties = _parties.size();
  // in-active parties's bit is set to true
  vector<bool> recvBits(nParties, true);

  // -- each party sends happiness to all others
  // -- set to unhappy if sees any unhappy
  _baParty.exchangeBit(happiness, recvBits);
  for (auto b : recvBits) {
    happiness &= b;
  }
  // -- all parties run consensus on happiness
  happiness =_baParty.consensus(happiness);
  return happiness;
}


// TODO: implement
template <class FieldType>
void LinearParty<FieldType>::
faultLocalizeConsistency(vector<int>& elimIds){
  cout << "fault localization for checking consistency not implemented yet"  << endl;
  assert(false);
  return;
}


// TODO implement
template <class FieldType>
void LinearParty<FieldType>::
faultLocalizeEvalSeg(vector<int>& elimIds){
  cout << "fault localization for eval seg not implemented yet"  << endl;
  assert(false);
  return;
}

// checks values spread by kingId are consistent
template <class FieldType>
bool LinearParty<FieldType>::
checkConsistency(int kingId, vector<FieldType> &vals,
                 vector<int>& elimIds) {
  bool happiness = true;
  int smallTp = (_nActiveParties - _bigT) / 2;
  int nPartiesInc = _parties.size() +1;
  int nVerifiers = _bigT + smallTp;

  if (vals.size() == 0) {
    return true;
  }

  // computational phase
  // -- all parties compute r^1 ... r^{T+t'} using HIM
  vector<FieldType> paddedVals = vals;
  paddedVals.resize(nPartiesInc, *(_field->GetZero()));
  vector<FieldType> expandedShares(nPartiesInc);
  _M.MatrixMult(paddedVals, expandedShares);
  expandedShares.resize(nVerifiers);

  // -- all parties send their values to corresponding dealer 
  vector<int> dealers;
  firstNActiveParties(nVerifiers, _activeMask, dealers);
  vector<FieldType> recvElms;
  gatherForDealers(expandedShares, recvElms, dealers);

  if (recvElms.size() > 0) {
    // I'm a verifier
    for (int i = 0; i < nPartiesInc - 1; i++) {
      happiness &= ( recvElms[i] == recvElms[i+1] );
    }
  }

  // fault detection phase
  bool success = faultDetection(happiness);

  if (!success) {
    faultLocalizeConsistency(elimIds);
  }

  return success;
}


template <class FieldType>
void LinearParty<FieldType>::
expandNShareForPartition(int partitionIdx,
                         vector<FieldType>& rShareNp,
                         vector<FieldType>& rShareT){
  int nPartiesInc = _parties.size() + 1;
  vector<FieldType> fTmp;
  vector<FieldType> alphaTmp(_smallT + 1);
  vector<FieldType> valTmp(_smallT + 1);

  int curInterpSize = _partition[ partitionIdx ].size();
  for (int i=0; i<curInterpSize; i++) {
    int id = _partition[0][i];
    alphaTmp[i] = _alpha[id];
    valTmp[i] = (_activeMask[id]) ? rShareNp[id] : _field->Random();
  }
  while (curInterpSize < _smallT + 1) {
    alphaTmp[ curInterpSize ] = // just an arbitrary value
      _field->GetElement( nPartiesInc + curInterpSize );
    valTmp[ curInterpSize ] = _field->Random();
    curInterpSize++;
  }
  _eccAlpha.interpolate(alphaTmp, valTmp, fTmp);
  rShareT.resize(nPartiesInc);
  for (int i = 0; i < nPartiesInc; i++) {
    rShareT[i] = _eccAlpha.evalPolynomial(_alpha[i], fTmp);
  }
  
  return;
}


// expand a whole (n'-1)-share into 3 r-shares
template <class FieldType>
void LinearParty<FieldType>::
expandNShare(vector<FieldType>& rShareNp,
             vector<FieldType>& rShareT1,
             vector<FieldType>& rShareT2,
             vector<FieldType>& rShareT3){

  expandNShareForPartition(0, rShareNp, rShareT1);
  expandNShareForPartition(1, rShareNp, rShareT2);
  expandNShareForPartition(2, rShareNp, rShareT3);
  
  return;
}

template <class FieldType>
void LinearParty<FieldType>::
extractNShare(vector<FieldType>& rShareT1,
              vector<FieldType>& rShareT2,
              vector<FieldType>& rShareT3,
              vector<FieldType>& rShareNp){

  int nPartiesInc = _parties.size() +1;
  rShareNp.resize();
  
  for (auto id : _partition[0]) {
    rShareNp[ id ] = rShareT1[id];
  }

  for (auto id : _partition[1]) {
    rShareNp[ id ] = rShareT2[id];
  }
  
  for (auto id : _partition[2]) {
    rShareNp[ id ] = rShareT3[id];
  }

  return;
}

template <class FieldType>
void LinearParty<FieldType>::
eliminateParties(const vector<int> &elimIds) {
  for (auto id : elimIds) {
    _activeMask[id] = false;
    _nActiveParties -= 1;
    _baParty.setPartyInactive(id);
    _eccAlpha.setPartyInactive(id);
    _eccBeta.setPartyInactive(id);
  }
  return;
}

// generate _bigT random 4-consistent tuples or return false
// Called QuadrupleShareRandom in the new paper
template <class FieldType>
bool LinearParty<FieldType>::
generate4Tuples(vector<vector<FieldType>> &rTuples) {
  bool happiness = true;
  int nPartiesInc = _parties.size() +1;
  FieldType r = _field->Random();
  rTuples.resize(4);

  // generate the 0_r_t share
  FieldType r0;
  happiness &= singleShareSecrete(r, _smallT, rTuples[0], r0);

  // -- locally generate the 1_r_t, 2_r_t, 3_r_t shares
  vector<FieldType> f;
  vector<FieldType> fVals(nPartiesInc);
  randomSecretePoly(r, _nActiveParties-1, f);
  for (int i = 0; i < nPartiesInc; i++) {
    fVals[i] = _eccAlpha.evalPolynomial(_alpha[i], fVals);
  }
  vector<FieldType> fValsT1, fValsT2, fValsT3;
  expandNShare(fVals, fValsT1, fValsT2, fValsT3);

  // -- share the 3 expanded vectors to all
  // -- apply HIM to received shares
  vector<FieldType> recvSharesT(nPartiesInc);
  vector< vector<FieldType> > transSharesT(3, vector<FieldType>(nPartiesInc));
  scatterForAll(fValsT1, recvSharesT);
  _M.MatrixMult(recvSharesT, transSharesT[0]);
  scatterForAll(fValsT2, recvSharesT);
  _M.MatrixMult(recvSharesT, transSharesT[1]);
  scatterForAll(fValsT3, recvSharesT);
  _M.MatrixMult(recvSharesT, transSharesT[2]);

  // -- every party send transformedShares of T+1 ... n' to
  // -- the corresponding parties
  int nVerifiers = _nActiveParties - _bigT;
  vector<int> dealers(nVerifiers);
  lastNActiveParties(nVerifiers, _activeMask, dealers);
  vector<FieldType> verifySharesT(nPartiesInc);
  vector< vector<FieldType> > toVerifySharesT(3);
  for (int i = 0; i < 3; i++) {
    lastNActiveValues(nVerifiers, _activeMask, transSharesT[0], verifySharesT);
    gatherForDealers(verifySharesT, toVerifySharesT[i], dealers);
  }

  if (toVerifySharesT[0].size() > 0) { // I'm a verifyer
    vector<FieldType> g;
    happiness &= checkPolynomial(toVerifySharesT[0], _smallT, g);
    happiness &= checkPolynomial(toVerifySharesT[1], _smallT, g);
    happiness &= checkPolynomial(toVerifySharesT[2], _smallT, g);
    
    vector<FieldType> rShareNp;
    extractNShare(toVerifySharesT[0],
                  toVerifySharesT[1],
                  toVerifySharesT[2], rShareNp);
    happiness &= checkPolynomial(rShareNp, _nActiveParties-1, g);
    happiness &= (_eccAlpha.evalPolynomial(*(_field->GetZero()), g) == r0);
  }

  // output the first T as results
  firstNActiveValues(_bigT, _activeMask, transSharesT[0], rTuples[1]);
  firstNActiveValues(_bigT, _activeMask, transSharesT[1], rTuples[2]);
  firstNActiveValues(_bigT, _activeMask, transSharesT[2], rTuples[3]);

  return happiness;
}


// TODO: implement
template <class FieldType>
bool LinearParty<FieldType>::
check4Consistency() {
  cout << "check4consistency not implemented yet"  << endl;
  assert(false);
  return true;
}

template <class FieldType>
bool LinearParty<FieldType>::
evalSeg(vector<TGate>& circuitSeg, vector<int>& elimIds){
  bool happiness = true;

  // -- generate enough multiplication tuples for current segment
  int multOffset = 0;
  vector< vector<FieldType> > a;
  vector< vector<FieldType> > b;
  vector< FieldType > c;
  // TODO, change to generateTuplesPE when implemented
  // and remove the assert().
  happiness &= generateTuples(a, b, c);
  assert(happiness);


  // if (!happiness) {
  //   return false;
  // }

  // -- For each addition gate: just add the two shares
  // -- For each mult gate: use a multiplication tuple (using lowest active id)
  // ---- Need to record the distribued d, and e values by kingId!
  // ---- Need to record the used xShare and yShare by mult gates!
  int kingId = 0;
  while(!_activeMask[kingId]){kingId++;}
  vector<FieldType> kingDs(_bigT, _zero);
  vector<FieldType> kingEs(_bigT, _zero);
  vector<FieldType> xShares(_bigT, _zero);
  vector<FieldType> yShares(_bigT, _zero);

  for (auto gate : circuitSeg) {
    FieldType xShare = _wireShares[gate.input1];
    FieldType yShare = _wireShares[gate.input2];
    switch (gate.gateType) {
    case SCALAR_ADD :
    case ADD :
      _wireShares[gate.output] = xShare + yShare;
      break;
    case SUB :
      _wireShares[gate.output] = xShare - yShare;
      break;
    case SCALAR :
    case MULT :
      xShares[multOffset] = xShare;
      yShares[multOffset] = yShare;
      _wireShares[gate.output] =
        evalMultGate(kingId, xShare, yShare,
                     a[0][multOffset], a[1][multOffset],
                     b[0][multOffset], b[1][multOffset], c[multOffset],
                     kingDs[multOffset], kingEs[multOffset]);
      multOffset++;
      break;
    default :
      cout << "Gate type Not Implemented." << endl;
      assert(false);
      break;
    }
  }

  if(multOffset == 0){
    return true;
  }
  
  // assert(multOffset == _bigT);
  kingDs.resize(multOffset);
  kingEs.resize(multOffset);
  // -- check the consistency of kingId
  happiness &= checkConsistency(kingId, kingDs, elimIds);
  if (!happiness) {
    return false;
  }
  happiness &= checkConsistency(kingId, kingEs, elimIds);
  if (!happiness) {
    return false;
  }

  // -- recompute all reconstructions (using recorded xShares and yShares)
  vector<FieldType> dShares(multOffset, _zero);
  vector<FieldType> eShares(multOffset, _zero);
  for (int i = 0; i < multOffset; i++) {
    dShares[i] = xShares[i] + a[0][i];
    eShares[i] = yShares[i] + b[0][i];
  }

  // -- if all check passes, return true.
  // -- otherwise, run faultLocalizationEvalSeg(). (= step 6 - 8 in paper)
  // -- TODO: maybe change here after implementing fault localization
  vector<FieldType> reconsDs(multOffset, _zero);
  happiness &= reconstructPublic(dShares, reconsDs, _smallT);
  // assert(happiness); // ^^  this is robust
  for (int i = 0; i < multOffset; i++) {
    if (reconsDs[i] != kingDs[i]) {
      faultLocalizeEvalSeg(elimIds);
      return false;
    }
  }
  vector<FieldType> reconsEs(multOffset, _zero);
  happiness &= reconstructPublic(eShares, reconsEs, _smallT);
  // assert(happiness); // ^^  this is robust
  for (int i = 0; i < multOffset; i++) {
    if (reconsEs[i] != kingEs[i]) {
      faultLocalizeEvalSeg(elimIds);
      return false;
    }
  }

  return happiness;
}

// fill a list of partyIds with inputs.
// also return the minimum number of input amoung them
static inline int
collectInputIds(vector<int> &inputSizes,
                vector<int> &inputParties,
                int maxBatch) {

  int nPartiesInc = inputSizes.size();
  inputParties.clear();
  
  // -- find first party with input
  int firstInputParty = 0;
  while (firstInputParty < nPartiesInc &&
         inputSizes[firstInputParty] == 0) {
    firstInputParty++;
  }
  if (firstInputParty == nPartiesInc) {
    return 0;
  }// else, found the first party with input

  
  // -- fill the list of Ids, and record minInputs.
  int minInputs = maxBatch;
  for (int i = firstInputParty; i < nPartiesInc; i++) {
    if (inputSizes[i] > 0) {
      inputParties.push_back(i);
      minInputs = min(minInputs, inputSizes[i]);
    }
  }
  
  // -- remove minInput from all
  for (int i = firstInputParty; i < nPartiesInc; i++) {
    if (inputSizes[i] > 0) {
      inputSizes[i] -= minInputs;
    }
  }
  return minInputs;
}

template <class FieldType>
static inline void
storeRecvDealerInputs(vector<FieldType>& totalInputs,
                      vector<FieldType>& dealerInputs,
                      const vector<int>& dealerIds) {
  int nDealers = dealerIds.size();
  for (int i = 0; i < nDealers; i++) {
    totalInputs[ dealerIds[i] ].insert( totalInputs[ dealerIds[i] ].end(),
                                        dealerInputs[i].begin(),
                                        dealerInputs[i].end() );
  }
  return;
}


// From BTH08 Appendex: batch input sharing (K dealers each T inputs)
// -- assuming everyone already have enough random sharings
// -- For K*T times: ReconsPriv() a random share to corresponding Dealer.
// -- Each Dealer computes T differences
// -- All dealers batch broadcast its T differences
// -- Each party locally compute input shares
template <class FieldType>
void LinearParty<FieldType>::
linearInput(vector< vector<FieldType> >& recvInputs,
            int nInputs,  vector<int>& dealerIds) {

  int nDealers = dealerIds.size();
  int nPartiesInc = _parties.size() + 1;
  int msgSize = _field->getElementSizeInBytes();
  recvInputs.resize(nDealers);
  sort(dealerIds.begin(), dealerIds.end());

  // -- store random shares in recvInputs first
  // -- later will add ``offsets'' to it
  for (int i = 0; i < nDealers; i++) {
    recvInputs[i].resize(nInputs);
    for (int j = 0; j < nInputs; j++) {
      recvInputs[i][j] = _singleRandomShares[_singleRandomSharesOffset++];
    }
  }

  // -- For K*T times: ReconsPriv() a random share to corresponding Dealer.
  // -- Each Dealer computes T differences
  vector<FieldType> sendOffsets(nInputs);
  for (int i = 0; i < nDealers; i++) {
    int dealerId = dealerIds[i];
    for (int j = 0; j < nInputs; j++) {
      FieldType result = 
        reconstructPrivate(recvInputs[i][j], _smallT, dealerId);

      if (dealerId == _myId) {
        sendOffsets[j] = _input[_inputOffset++] - result;
      }
    }
  }

  // -- All dealers batch broadcast its T differences
  vector< vector<FieldType> > recvOffsets(nPartiesInc);
  robustBatchBroadcast(sendOffsets, recvOffsets, nInputs, dealerIds);

  // -- Each party locally compute input shares
  for (int i = 0; i < nDealers; i++) {
    for (int j = 0; j < nInputs; j++) {
      recvInputs[i][j] += recvOffsets[i][j];
    }
  }

  return;
}

// runs linearInput() repeatedly until all inputs are shared.
// TODO: also plug in random shares for random gate
template <class FieldType>
void LinearParty<FieldType>::
InputPhase(){

  _wireShares.resize( _circuit.getNrOfGates(), _zero );
  int nPartiesInc = _parties.size() + 1;
  
  // batch share inputs until all inputs are shared
  vector< vector<FieldType> > receivedInputShares(nPartiesInc);
  vector<int> inputParties;
  vector<int> inputsToSend = _inputSizes;
  int minInputs = collectInputIds(inputsToSend, inputParties, _bigT);
  while(minInputs > 0){
    vector< vector<FieldType> > batchRecvInput;
    linearInput(batchRecvInput, minInputs, inputParties);
    storeRecvDealerInputs(receivedInputShares, batchRecvInput, inputParties);
    minInputs = collectInputIds(inputsToSend, inputParties, _bigT);
  }

  cout << "-------- finished input sharing" << endl;

  // // batch reconstruct inputs until all inputs are reconstructed
  // vector< vector<FieldType> > receivedInputsClear(nPartiesInc);
  // vector<FieldType> reconsShares;
  // vector<FieldType> reconsResults;
  // for (int i = 0; i < nPartiesInc; i++) {
  //   int nInputs = receivedInputShares[i].size();
  //   receivedInputsClear[i].resize( nInputs );
  //   if (nInputs == 0) {
  //     continue;
  //   }

  //   for (int j = 0; j < nInputs; j += _bigT) {
  //     int nShares = (nInputs - j) < _bigT ? (nInputs - j) : _bigT;
            
  //     reconsShares.resize(nShares);
  //     for (int k = 0; k < nShares; k++) {
  //       reconsShares[k] = receivedInputShares[i][j+k];
  //     }

  //     reconstructPublic(reconsShares, reconsResults, _smallT);
      
  //     for (int k = 0; k < nShares; k++) {
  //       receivedInputsClear[i][j+k] = reconsResults[k];
  //     }
  //   }
  // }

  // cout << "-------- finished input reconstruction" << endl;

  // plug in shares of  input values
  int nInputGates = _circuit.getNrOfInputGates();
  const vector<TGate> gates = _circuit.getGates();
  vector<int> offsets(nPartiesInc);

  int i=0;
  while (i < nInputGates) {
    if (gates[i].gateType == INPUT) {
      int curParty = gates[i].party;
      int curOffset = offsets[ curParty ]++;
      _wireShares[ gates[i].output ] = receivedInputShares[curParty][curOffset];
      i++;
    }
  }
  cout << "---- finished Input Phase ----" << endl;
  return;
}


// evaluate circuits
// -- repeatedly run evalSeg until all segments are done.
// -- eliminate players if evalSeg fails, and re-do.
template <class FieldType>
void LinearParty<FieldType>::
EvalPhase(){  
  const vector<TGate> gates = _circuit.getGates();
  // TODO: should also subtract random gates here?
  int nTotalGates = _circuit.getNrOfGates();
  int nGatesLeft = nTotalGates -
    _circuit.getNrOfInputGates() - _circuit.getNrOfOutputGates();
  int gateIdx = 0;
  int nMults = 0;
  
  vector<TGate> seg;
  while (nGatesLeft > 0) {
    // build a segment
    seg.clear();
    nMults = 0;
    while (nMults < _bigT && gateIdx < nTotalGates) {
      auto g = gates[gateIdx++];
      if (g.gateType != INPUT && g.gateType != OUTPUT) {
        seg.push_back( g );
        nMults += (g.gateType == MULT || g.gateType == SCALAR);
        nGatesLeft--;
      }
    }
    // cout << "evaluating " << gateIdx << "/" << nTotalGates
    //      << " with " << nMults << " mult gates "<< endl;
    // evaluate a segment (repeatedly until success)
    vector<int> elimIds; // the first time, just eliminate nothing.
    bool success = false;
    while (!success) {
      eliminateParties(elimIds);
      success = evalSeg(seg, elimIds);
    }
  }
  
  cout << "---- finished Eval Phase ----" << endl;
}


// reconstruct outputs
template <class FieldType>
void LinearParty<FieldType>::
OutputPhase(){
  int nOutputGates = _circuit.getNrOfOutputGates();
  const vector<TGate> gates = _circuit.getGates();

  cout << "---- there are " << nOutputGates
       << " output gates" << endl;

  for (auto g : gates) {
    if (g.gateType == OUTPUT) {
      // for each output gate, consume a zero-t-share
      int rootId = g.party;
      FieldType outputShare = _wireShares[g.input1];
      outputShare += _singleZeroShares[_singleZeroSharesOffset++];
      FieldType output = reconstructPrivate(outputShare, _smallT, rootId);
      if (rootId == _myId) {
        cout << _myId << " : wire " << g.input1 << " is "
             << output << endl;
      }
    }
  }
  
  cout << "---- finished Output Phase ----" << endl;
  return;
}


// runs preparation phase:
// -- generate random shares for input and random gates
// -- generate (random) zero shares for output gates
// -- TODO: fault detection, localization, and re-run
template <class FieldType>
void LinearParty<FieldType>::
runOffline(){
  // -- generate random t-shares for input and random gates
  int nSingleShares =
    _circuit.getNrOfInputGates() + _circuit.getNrOfRandomGates();
  int nBatches = (nSingleShares + _bigT-1) / _bigT;

  cout << "running singleShareRandom in "
       << nBatches << " batches" << endl;
  
  _singleRandomShares.resize(nBatches * _bigT);
  vector<FieldType> batchResult;
  bool happiness = true;
  for (int i = 0; i < nBatches; i++) {
    FieldType tmp;
    happiness &= singleShareRandom(_smallT, batchResult, tmp);
    for (int j = 0; j < _bigT; j++) {
      _singleRandomShares[ i*_bigT + j ] = batchResult[j];
    }
  }

  // -- generate (random) zero t-shares for output gates
  int nSingleZeroShares = _circuit.getNrOfOutputGates();
  nBatches = (nSingleZeroShares + _bigT-1) / _bigT;

  cout << "running singleShareZero in "
       << nBatches << " batches" << endl;
  
  _singleZeroShares.resize(nBatches * _bigT);
  happiness = true;
  for (int i = 0; i < nBatches; i++) {
    happiness &= singleShareZero(_smallT, batchResult);
    for (int j = 0; j < _bigT; j++) {
      _singleZeroShares[ i*_bigT + j ] = batchResult[j];
    }
  }
  
  cout << "---- finished offline phase ----" << endl;
  return;
}


// runs =main()= protocol
// -- InputPhase(): share inputs
// -- EvalPhase(): evaluate circuit
// -- OutputPhase(): reconstruct outputs
template <class FieldType>
void LinearParty<FieldType>::
runOnline(){
  InputPhase();
  EvalPhase();
  OutputPhase();
  return;
}

#endif /* LINEARPARTY_H_ */

