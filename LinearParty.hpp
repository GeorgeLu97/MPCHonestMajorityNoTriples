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
// - circuit reading/ calculation
// - BA broadcast and (batch) Input sharing
// - 4 consistency functionality
//
// - reconstruction
// - fault localization
// - player elimination
//
// - a running protocol (with overall procedure =eval()=)
// - performance testing
// - (trivial) malicious party
// - optimizations
//

// TODO:
// 1. use broadcast functionality from BAParty
// 2. 4 consistency functionality
// 3. reconstruction
// 4. fault localization
template <class FieldType>
class LinearParty : public Protocol, public HonestMajority, MultiParty{
  
private:
  
  // ======== protocol properties ========
  int _myId;
  int _nActiveParties;
  int _nThread;
  ArithmeticCircuit _circuit;
  TemplateField<FieldType>*_field;

  // manage active/inactive parties by a mask
  vector< shared_ptr<ProtocolPartyData> > _parties;
  vector<bool> _activeMask;
  
  
  // only store long value (not FieldType) to save memory
  vector<int> _inputSizes;
  vector<long> _input;
  vector<long> _output;


  // ======== private functionalities ========
  
  // ---- helper funtions ----
  inline string getArg(string argKey)
  { return this->getParser().getValueByKey(this->arguments, argKey); }
  void makeField();
  void makeParties();

  void encodeFieldElts(vector<FieldType>& input, vector<byte>& output);
  void decodeFieldElts(vector<byte>& input, vector<FieldType>& output);  

  void broadcastWorker(vector<byte>& msg, int threadId, int nThreads);
  void broadcastElements(vector<FieldType>& buffer,
                         int nElements, int rootId);

  // ---- main subprotocols ----
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
  vector<FieldType> _wireValues;
    
};

static inline void 
readInputFromFile(vector<long>& input, string fileName){

  ifstream inputFile(fileName);
  int inputSize = input.size();

  // read available inputs to fill _input
  for(int i=0; i< inputSize && !(inputFile.eof()); i++){
    inputFile >> input[i];
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
  
  // build _circuit
  _circuit.readCircuit(getArg("circuitFile"));
  _circuit.reArrangeCircuit();
  _input.resize(_circuit.getNrOfInputGates());
  _output.resize(_circuit.getNrOfOutputGates());
  
  // build _inputSizes
  fillInputSizes(_myId, _parties.size() + 1, _circuit, _inputSizes);
  readInputFromFile(_input, getArg("inputFile"));

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
  
}

template <class FieldType>
void LinearParty<FieldType>::
broadcastWorker(vector<byte>& msg,
                int threadId, int nThreads) {
  int nParties = _parties.size();
  for (int i = threadId; i < nParties; i += nThreads) {
    if (!_activeMask[i]) {
      continue;
    }

    _parties[i]->getChannel()->write(msg.data(), msg.size());
  }
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
decodeFieldElts(vector<byte>& input, vector<FieldType>& output){
  int fieldByteSize = _field->getElementSizeInBytes();
  // if(input.size() % fieldByteSize) { cout << "input byte misalignment" << endl; }
  output.resize(input.size() / fieldByteSize);

  for(int j = 0; j < output.size(); j++) {
    output[j] = _field->bytesToElement(input.data() + fieldByteSize * j);
  }
}


// broadcast from a root to all
// a specified number of field elements
// -- TODO: replace with BAParty functionality
// -- TODO: currently assumes a broadcast channel.
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
    
    vector<thread> threads(_nThread);
    for (int i = 0; i < _nThread; i++) {
      threads[i] = thread(&LinearParty::broadcastWorker, this,
                          ref(msg), i, _nThread);
    }
    for(int i=0; i<_nThread; i++){
      threads[i].join();
    }
    
  } else {
    // I'm a receiver
    msg.resize(msgSize);
    for (int i = 0; i < nParties; i++) {
      if (!_activeMask[i] ||
          _parties[i]->getID() != rootId) {
        continue;
      }
      _parties[i]->getChannel()->read(msg.data(), msg.size());
    }
    decodeFieldElts(msg, buffer);
  }
  return;
}



// share inputs
// -- TODO: currently just a stub: share inputs in the clear
template <class FieldType>
void LinearParty<FieldType>::
InputPhase(){
  _wireValues.resize( _circuit.getNrOfGates() );

  // broadcast for all parties
  int nPartiesInc = _parties.size() + 1;
  vector< vector<FieldType> > receivedInputs(nPartiesInc);
  for(int i=0; i<nPartiesInc; i++){
    int curInputSize = _inputSizes[i];
    if (_myId == i) {
      // I'm the root
      receivedInputs[i].resize(curInputSize);
      for (int j = 0; j < curInputSize; j++) {
        receivedInputs[j] = FieldType(_input[j]);
      }
    }
    broadcastElements(receivedInputs[i], curInputSize, i);
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
    }
  }

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
      _wireValues[g.output] = _wireValues[g.input1] + _wireValues[g.input2];
      break;
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
}

// create random shares
// -- TODO: currently a stub: do nothing.
template <class FieldType>
void LinearParty<FieldType>::
RandomPhase(){
  return;
}



// reconstruct outputs
// -- TODO: currently just a stub: do nothing.
template <class FieldType>
void LinearParty<FieldType>::
OutputPhase(){
  
  int nOutputGates = _circuit.getNrOfOutputGates();
  const vector<TGate> gates = _circuit.getGates();

  int i=0;
  while (i < nOutputGates) {
    if (gates[i].gateType == OUTPUT && gates[i].party == _myId) {
      cout << _myId << " : output " << gates[i].output
           << " is " << _wireValues[ gates[i].output ];
    }
  }
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
