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
// - 4 consistency functionality
// - fault localization
// - player elimination
// 
// - performance testing
// - (trivial) malicious party
// - optimizations
//

// TODO:
// 1.4 check that all inactive elements are 0.

// 2. 4 consistency functionality
// 4. fault localization

// Questions:
// 1. Error checking during protocol (e.g. timeout communication?)
// 2. How to agree on a constant HIM matrix (currently just 1 ~ 2n)?
// 3. updating HIM or not?

template <class FieldType>
class LinearParty : public Protocol, public HonestMajority, MultiParty{
  
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

  vector<int> _inputSizes;
  vector<long> _output;

  int _inputOffset = 0;
  vector<FieldType> _input;
  int _singleRandomSharesOffset = 0;
  vector<FieldType> _singleRandomShares;
  int _singleZeroSharesOffset = 0;
  vector<FieldType> _singleZeroShares;
  int _tripleRandomSharesOffset = 0;
  vector< vector<FieldType> >_tripleRandomShares;


  // ======== private functionalities ========
  
  // ---- helper funtions ----
  inline string getArg(string argKey)
  { return this->getParser().getValueByKey(this->arguments, argKey); }
  void encodeFieldElts(vector<FieldType>& input, vector<byte>& output);
  void encodeFieldElt(FieldType& input, vector<byte>& output);
  void decodeFieldElts(vector<byte>& input, vector<FieldType>& output);
  void decodeFieldElt(vector<byte>& input, FieldType& output);
  void readInputFromFile(string fileName);
  void randomSecretePoly(FieldType& s, int degree,
                         vector<FieldType>& polynomial);
  void firstNActiveParties(int N, vector<int>& activeIds);
  void lastNActiveParties(int N, vector<int>& activeIds);

  // ---- helper round functions ----
  void scatterForDealers(vector<FieldType>& sendElms,
                         vector<FieldType>& recvElms,
                         vector<int>& dealerIds);
  void gatherForDealers(vector<FieldType>& sendElms,
                        vector<FieldType>& recvElms,
                        vector<int>& dealerIds);
  void scatterForAll(vector<FieldType>& sendElms,
                     vector<FieldType>& recvElms);

  // ---- initialization funtions ----
  void makeField();
  void makeParties();
  void makeAlpha();
  void makeBeta();

  // ---- Batch broadcasting (BTH08 Appendix) ----
  // expand polynomial to n evaluations. Encode each into messages
  // also return the evaluation on my own alpha.
  void prepareExpandedElms(vector<FieldType>& polynomial,
                           vector<FieldType>& expandedElms);
  // Broadcast (each dealer (k intotal) spreads T values)
  void robustBatchBroadcast(vector<FieldType>& sendElms,
                            vector< vector<FieldType> >& recvElms,
                            int nElements,  vector<int>& dealerIds);
  // robust reconstruct nInputs # of input.
  // -- only works when nInputs <= _bigT
  void linearInput(vector< vector<FieldType> >& recvInputs,
                   int nInputs,  vector<int>& dealerIds);
  

  // ---- Random shares and Reconstructions ----
  // create T random shares among parties. return false if fails
  // -- For a verifier party in the protocol,
  // -- ``checkVal'' is the reconstructed polinomial evaluated at 0.
  // -- For others, just 0.
  // Note: non-robust.
  bool singleShareRandom(int d, vector<FieldType>& shares,
                         FieldType& checkVal);
  // similarly create T random multiple-shares among parties
  bool multipleShareRandom(const vector<int> degrees,
                           vector< vector<FieldType> >& shares);
  // similarly create T random shares of 0 among parties
  bool singleShareZero(int d, vector<FieldType> shares);

  // reconstruct an element towards root. return 0 if fails
  // Note: robust if degree <= T. 
  FieldType reconstructPrivate(FieldType& share, int degree, int root);
  // reconstruct (at most) T elements towards all. return false if fails.
  // Note: robust if degree <= T.
  bool reconstructPublic(vector<FieldType>& shares,
                         vector<FieldType>& results, int degree);

  // ---- main subprotocols ----
  // From BTH08 Appendix: batch input sharing (K dealers each T inputs)
  void InputPhase();
  // From new paper =eval()=
  void EvalPhase();
  // From new paper =output()=
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
  void broadcastElements(vector<FieldType>& buffer,
                         int nElements, int rootId);

  vector<FieldType> _wireValues;
    
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

  // initialize subprotocol classes
  _M.allocate(_nActiveParties, _nActiveParties, _field);
  _M.InitHIMByVectors(_alpha, _beta);
  _baParty.setParties(_parties, _myId);
  _baParty.setNumThreads(_nThread);
  _baParty.setAlphaBeta(_alpha, _beta);
  _baParty.setHIM(_M);
  _eccAlpha.setAlpha(_alpha);
  _eccAlpha.setField(_field);
  _eccBeta.setAlpha(_beta);
  _eccBeta.setField(_field);
      
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


// broadcast from a root to all
// a specified number of field elements
template <class FieldType>
void LinearParty<FieldType>::
broadcastElements(vector<FieldType>& buffer,
                  int nElements, int rootId){

  vector<byte> msg;
  int nParties = _parties.size();
  int msgSize = _field->getElementSizeInBytes() * nElements;

  if (rootId == _myId) {
    // I'm the sender
    encodeFieldElts(buffer, msg);
    msg.resize(msgSize);
    _baParty.broadcastMsg(msg, msgSize, rootId);
  } else {
    // I'm a receiver
    _baParty.broadcastMsg(msg, msgSize, rootId);
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
gatherForDealers(vector<FieldType>& sendElms,
                 vector<FieldType>& recvElms,
                 vector<int>& dealerIds) {
  int nDealers = dealerIds.size();
  int nParties = _parties.size();
  int nPartiesInc = _parties.size() +1;

  recvElms.resize(nPartiesInc);
  sendElms.resize(nDealers);
  
  int msgSize = _field->getElementSizeInBytes();
  vector<byte> sendMsg(msgSize);
  vector< vector<byte> > recvMsgs(nParties, vector<byte>(msgSize));
  for(int i=0; i<nDealers; i++){
    int dealerId = dealerIds[i];

    if ( dealerId == _myId) {
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
void LinearParty<FieldType>::
lastNActiveParties(int lastN, vector<int>& activeIds) {
  int nPartiesInc = _activeMask.size();
  int lastId = nPartiesInc-1;
  activeIds.resize(lastN);
  while (lastN > 0 && lastId >= 0) {
    if (_activeMask[lastId]) {
      activeIds[lastN-1] = lastId;
      lastN--;
    }
    lastId--;
  }
  return;
}

template <class FieldType>
void LinearParty<FieldType>::
firstNActiveParties(int firstN, vector<int>& activeIds) {
  int nPartiesInc = _activeMask.size();
  int firstId = 0;
  int idx = 0;
  activeIds.resize(firstN);
  while (idx < firstN) {
    if (_activeMask[firstId]) {
      activeIds[idx] = firstId;
      idx++;
    }
    firstId++;
  }
  return;
}


// create T random shares among parties. return false if fails
// Note: non-robust.
// -- every party choose a random secrete, and share with all others.
// -- every party applies HIM on received shares
// -- every party send transformedShares of T+1 ... n' to crrspnd. parties
// -- the (n'-T) verifiers verify (n'-T) received shares
template <class FieldType>
bool LinearParty<FieldType>::
singleShareRandom(int d, vector<FieldType>& shares, FieldType& checkVal) {
  bool happiness = true;
  int nPartiesInc = _parties.size() + 1;

  // -- every party choose a random secrete
  FieldType r = _field->Random();
  vector<FieldType> f;
  vector<FieldType> sendShares(nPartiesInc);
  randomSecretePoly(r, d, f);
  for (int i = 0; i < nPartiesInc; i++) {
    sendShares[i] = _eccAlpha.evalPolynomial(_alpha[i], f);
  }
  // -- and share with all others.
  vector<FieldType> recvShares(nPartiesInc);
  scatterForAll(sendShares, recvShares);

  // -- every party applies HIM on received shares
  vector<FieldType> transformedShares(nPartiesInc);
  _M.MatrixMult(recvShares, transformedShares);

  // -- every party send transformedShares of T+1 ... n' to
  // -- the corresponding parties
  int nVerifiers = _nActiveParties - _bigT;
  vector<int> dealers(nVerifiers);
  lastNActiveParties(nVerifiers, dealers);
  vector<FieldType> verifyShares(nVerifiers);
  vector<FieldType> toVerifyShares;
  for (int i = 0; i < nVerifiers; i++) {
    verifyShares[i] = transformedShares[ dealers[i] ];
  }
  gatherForDealers(verifyShares, toVerifyShares, dealers);

  // -- the (n'-T) verifiers verify (n'-T) received shares
  vector<FieldType> g;
  happiness &= _eccAlpha.reconstruct(toVerifyShares, d, g);
  assert(happiness);
  
  shares.resize(_bigT);
  vector<int> activeSharers(_bigT);
  firstNActiveParties(_bigT, activeSharers);
  for (int i = 0; i < _bigT; i++) {
    shares[i] = transformedShares[ activeSharers[i] ];
  }

  return happiness;
}

// similarly create T random multiple-shares among parties
template <class FieldType>
bool LinearParty<FieldType>::
multipleShareRandom(const vector<int> degrees,
                    vector< vector<FieldType> > &shares) {
  bool happiness = true;
  int nShares = degrees.size();
  int nPartiesInc = _parties.size() + 1;

  return happiness;
}

// similarly create T random shares of 0 among parties
// -- TODO implement
template <class FieldType>
bool LinearParty<FieldType>::
singleShareZero(int d, vector<FieldType> shares) {
  return true;
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
// -- TODO currently reconstructPublic all inputs after sharing
template <class FieldType>
void LinearParty<FieldType>::
InputPhase(){

  _wireValues.resize( _circuit.getNrOfGates() );
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
  
  // batch reconstruct inputs until all inputs are reconstructed
  vector< vector<FieldType> > receivedInputsClear(nPartiesInc);
  vector<FieldType> reconsShares;
  vector<FieldType> reconsResults;
  for (int i = 0; i < nPartiesInc; i++) {
    int nInputs = receivedInputShares[i].size();
    receivedInputsClear[i].resize( nInputs );
    if (nInputs == 0) {
      continue;
    }

    for (int j = 0; j < nInputs; j += _bigT) {
      int nShares = (nInputs - j) < _bigT ? (nInputs - j) : _bigT;
            
      reconsShares.resize(nShares);
      for (int k = 0; k < nShares; k++) {
        reconsShares[k] = receivedInputShares[i][j+k];
      }

      reconstructPublic(reconsShares, reconsResults, _smallT);
      
      for (int k = 0; k < nShares; k++) {
        receivedInputsClear[i][j+k] = reconsResults[k];
      }
    }
  }

  cout << "-------- finished input reconstruction" << endl;

  // plug in (reconstructed) clear input values
  int nInputGates = _circuit.getNrOfInputGates();
  const vector<TGate> gates = _circuit.getGates();
  vector<int> offsets(nPartiesInc);

  int i=0;
  while (i < nInputGates) {
    if (gates[i].gateType == INPUT) {
      int curParty = gates[i].party;
      int curOffset = offsets[ curParty ]++;
      _wireValues[ gates[i].output ] = receivedInputsClear[curParty][curOffset];
      i++;
    }
  }
  cout << "---- finished Input Phase ----" << endl;
  return;
}

// evaluate circuits
// -- TODO: currently just a stub: compute gates in clear
template <class FieldType>
void LinearParty<FieldType>::
EvalPhase(){  
  const vector<TGate> gates = _circuit.getGates();
  for (auto g : gates) {
    switch (g.gateType) {
    case ADD :
    case SCALAR_ADD :
      _wireValues[g.output] = _wireValues[g.input1] + _wireValues[g.input2];
      break;
    case SCALAR :
    case MULT :
      _wireValues[g.output] = _wireValues[g.input1] * _wireValues[g.input2];
      break;
    case SUB :
      _wireValues[g.output] = _wireValues[g.input1] - _wireValues[g.input2];
      break;
    default :
      break;
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
    if ((g.output == 8002 ||
         g.output == 8001 ||
         g.output == 5598 ||
         g.output == 5599 ||
         g.output == 7199 ||
         g.output == 7200)
        && _myId == 0) {
      cout << _myId << " : wire " << g.output
           << " is " << _wireValues[ g.output ] << endl;
      
    }

    if (g.gateType == OUTPUT && g.party == _myId) {
      cout << _myId << " : output " << g.input1
           << " is " << _wireValues[ g.input1 ] << endl;
    }
  }

  cout << "---- finished Output Phase ----" << endl;
  return;
}


// runs preparation phase:
// -- generate random shares for input and random gates
// -- TODO: generating shares random triples
template <class FieldType>
void LinearParty<FieldType>::
runOffline(){
  // -- generate random shares for input and random gates
  int nSingleShares =
    _circuit.getNrOfInputGates() + _circuit.getNrOfRandomGates();
  int nBatches = nSingleShares / _bigT;

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

