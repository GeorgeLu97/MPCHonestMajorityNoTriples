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
//
// - 4 consistency functionality
//
// - reconstruction
// - fault localization
// - player elimination
//
// - performance testing
// - (trivial) malicious party
// - optimizations
//

// TODO:
// 1. change all FieldType(n) to _field->getElement(n);
// -- GF2E will not properly initialize!!
// -- also do this for ECC.h and BAParty.h
// -- need to add _field member for them too.
// 2. 4 consistency functionality
// 3. reconstruction
// 4. fault localization

// Questions:
// 1. Error checking during protocol (e.g. timeout communication?)
// 2. How to agree on a constant HIM matrix?

template <class FieldType>
class LinearParty : public Protocol, public HonestMajority, MultiParty{
  
private:
  
  // ======== protocol properties ========
  int _myId;
  int _nActiveParties;
  int _nThread;
  ArithmeticCircuit _circuit;
  TemplateField<FieldType>*_field;
  
  // _alpha[i] is the field element of the party with id = i.
  // Note: every party has the same _alpha.
  //  ..   Same with _inputSizes
  vector<FieldType> _alpha;
  vector<int> _inputSizes;

  // subprotocol classes
  BAParty<FieldType> _baParty;
  ECC<FieldType> _ecc;

  // manage active/inactive parties by a mask
  vector< shared_ptr<ProtocolPartyData> > _parties;
  vector<bool> _activeMask;
  
  // only store long value (not FieldType) to save memory
  vector<FieldType> _input;
  vector<long> _output;


  // ======== private functionalities ========
  
  // ---- helper funtions ----
  inline string getArg(string argKey)
  { return this->getParser().getValueByKey(this->arguments, argKey); }
  void encodeFieldElts(vector<FieldType>& input, vector<byte>& output);
  void encodeFieldElt(FieldType& input, vector<byte>& output);
  void decodeFieldElts(vector<byte>& input, vector<FieldType>& output);
  void decodeFieldElt(vector<byte>& input, FieldType& output);
  void readInputFromFile(string fileName);

  // initialization funtions
  void makeField();
  void makeParties();
  void makeAlpha();

  // expand polynomial to n evaluations. Encode each into messages
  // also return the evaluation on my own alpha.
  FieldType prepareExpandedMsgs(vector<FieldType>& polynomial,
                                vector< vector<byte> >& expandedMsgs);
  // Broadcast (each dealer (k intotal) spreads T values)
  // -- Every dealer expand T values into n by interpolating
  // -- Every dealer distribute 1 value to each party
  // -- Run BroadcastForP() for k values (1 from each dealer)
  // -- Every player reconstruct T values for each dealer w/ ECC
  void robustBatchBroadcast(vector<FieldType>& sendElms,
                            vector< vector<FieldType> >& recvElms,
                            int nElements,  vector<int>& dealerIds);

  // ---- main subprotocols ----

  // From BTH08 Appendex: batch input sharing (K dealers each T inputs)
  void InputPhase();
  void RandomPhase();
  void EvalPhase();
  void OutputPhase();
  
public:

  LinearParty(int argc, char* argv[]);
  ~LinearParty();

  // ======== inherited functionalites ========
  bool hasOffline() override {return true;}
  bool hasOnline() override {return true;}
  void run() override;
  void runOffline() override;
  void runOnline() override;

  // ======== testing functionalities ========
  void broadcastElements(vector<FieldType>& buffer,
                         int nElements, int rootId);

  vector<FieldType> _wireValues;
    
};

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
  _activeMask.resize(_parties.size(), true);
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

template <class FieldType>
LinearParty<FieldType>::
LinearParty(int argc, char* argv[])
  : Protocol("MPCLinearComm", argc, argv)
{
  // fill basic class properties
  _myId = stoi(getArg("partyID"));
  _nActiveParties = stoi(getArg("partiesNumber"));
  _nThread = stoi(getArg("numThreads"));
  makeField();
  makeParties();
  makeAlpha();

  // initialize subprotocol classes
  _baParty.setParties(_parties, _myId);
  _baParty.setNumThreads(_nThread);
  _ecc.setAlpha(_alpha);
  
  
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


// runs preparation phase:
// -- TODO: generating random shares of triples
template <class FieldType>
void LinearParty<FieldType>::
runOffline(){
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
FieldType LinearParty<FieldType>::
prepareExpandedMsgs(vector<FieldType>& polynomial,
                    vector< vector<byte> >& expandedMsgs){

  int nParties = _parties.size();
  expandedMsgs.resize(nParties);
  for (int j = 0; j < nParties; j++) {
    FieldType currentPartyAlpha = _alpha[ _parties[j]->getID() ];
    FieldType currentEval = _ecc.evalPolynomial(currentPartyAlpha, polynomial);
    encodeFieldElt( currentEval, expandedMsgs[j]);
  }

  FieldType myAlpha = _alpha[ _myId ];
  return _ecc.evalPolynomial(myAlpha, polynomial);
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
  sort(dealerIds.begin(), dealerIds.end());

  // get one element from each dealer Id
  vector<FieldType> recvExpandedElms(nDealers);
  for(int i=0; i<nDealers; i++){
    vector< vector<byte> > sendMsgs;
    vector<byte> recvMsg;
    if (dealerIds[i] == _myId) {
      // prepare messages and scatter
      sendElms.resize(nElements);
      // store my own share, and build msgs to scatter
      recvExpandedElms[i] =  prepareExpandedMsgs(sendElms, sendMsgs);
      _baParty.scatterMsg(sendMsgs, recvMsg, elmSize, dealerIds[i]);
      
    } else {
      // receive message and store
      _baParty.scatterMsg(sendMsgs, recvMsg, elmSize, dealerIds[i]);
      decodeFieldElt(recvMsg, recvExpandedElms[i]);
    }
  }

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
    _ecc.reconstruct(allExpandedElms[i], nElements-1, recvElms[i]);
  }

  return;
}

// fill a list of partyIds with inputs.
// also return the minimum number of input amoung them
static inline int
collectInputIds(vector<int> &inputSizes,
                vector<int> &inputParties,
                int maxBatch) {

  int nPartiesInc = inputSizes.size();
  
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
InputPhase(){
    _wireValues.resize( _circuit.getNrOfGates() );
  int nPartiesInc = _parties.size() + 1;
  int smallT = _nActiveParties / 3;
  int bigT = _nActiveParties - 2*smallT;
  cout << "bigT is " << bigT << endl;
  
  // batch share inputs until all inputs are shared
  vector< vector<FieldType> > receivedInputs(nPartiesInc);
  vector<int> inputParties;
  vector<int> inputsToSend = _inputSizes;
  int minInputs = collectInputIds(inputsToSend, inputParties, bigT);
  auto inputHead = _input.begin();
  while(minInputs > 0){
    cout << "sending " << minInputs << " inputs" << endl;
    
    vector<FieldType> batchInput(inputHead, inputHead+minInputs);
    inputHead += minInputs;
    vector< vector<FieldType> > batchRecvInput;
    robustBatchBroadcast(batchInput, batchRecvInput, minInputs, inputParties);
    storeRecvDealerInputs(receivedInputs, batchRecvInput, inputParties);
    minInputs = collectInputIds(inputsToSend, inputParties, bigT);
  }


  // plug in input values
  int nInputGates = _circuit.getNrOfInputGates();
  const vector<TGate> gates = _circuit.getGates();
  vector<int> offsets(nPartiesInc);

  int i=0;
  while (i < nInputGates) {
    if (gates[i].gateType == INPUT) {
      int curParty = gates[i].party;
      int curOffset = offsets[ curParty ]++;
      _wireValues[ gates[i].output ] = receivedInputs[curParty][curOffset];
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

// create random shares
// -- TODO: currently a stub: do nothing.
template <class FieldType>
void LinearParty<FieldType>::
RandomPhase(){
  return;
}



// reconstruct outputs
template <class FieldType>
void LinearParty<FieldType>::
OutputPhase(){
  int nOutputGates = _circuit.getNrOfOutputGates();
  const vector<TGate> gates = _circuit.getGates();

  cout << "==== there are " << nOutputGates
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


// runs =main()= protocol
// -- TODO: 
// -- InputPhase(): share inputs
// -- RandomPhase(): share random gates
// -- EvalPhase(): evaluate circuit
// -- OutputPhase(): reconstruct outputs
template <class FieldType>
void LinearParty<FieldType>::
runOnline(){

  InputPhase();
  RandomPhase();
  EvalPhase();
  OutputPhase();
  return;
}




#endif /* LINEARPARTY_H_ */
