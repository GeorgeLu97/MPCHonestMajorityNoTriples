#ifndef LINEARPARTY_H_
#define LINEARPARTY_H_

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
// - circuit reading/ calculation/ writing 
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
template <class FieldType>
class LinearParty : public Protocol, public HonestMajority, MultiParty{
  
private:
  
  // -------- protocol properties --------
  int _myId;
  int _nParties;
  ArithmeticCircuit _circuit;
  TemplateField<FieldType>*_field;
  
  
  // only store long value (not FieldType) to save memory
  vector<int> _inputSizes;
  vector<long> _input;
  vector<long> _output;


  // -------- private functionalities --------
  static inline string getArg(string argKey)
  { return this->getParser().getValueByKey(this->arguments, argKey); }
  void makeField();
public:

  LinearParty(int argc, char* argv[]);
  ~LinearParty();
  
};


static inline void 
readInputFromFile(vector<long>& input, string fileName){

  ifstream inputFile(fileName);
  int inputSize = input.size();

  // read available inputs to fill _input
  for(int i=0; i< && !(inputFile.eof()); i++){
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
  for(int i=0; i<nInputGates; i++){
    if(gates[i].gateType == INPUT){
      sizes[gates[i].party]++;
    }
  }
  
  return;
}

template<class FieldType>
void LinearParty<class FieldType>::
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

template <class FieldType>
LinearParty<FieldType>::
LinearParty(int argc, char* argv[])
  : Protocol("MPCLinearComm", argc, argv)
{
  // fill basic class properties
  _myId = stoi(getArg("partyID"));
  _nParties = stoi(getArg("partiesNumber"));
  makeField();
  
  // build _circuit
  _circuit.readCircuit(getArg("circuitFile"));
  _circuit.reArrangeCircuit();
  _input.resize(_circuit.getNrOfInputGates());
  _output.resize(_circuit.getNrOfOutputGates());
  
  // build _inputSizes
  fillInputSizes(_myId, _nParties, _circuit, _inputSizes);
  readInputFromFile(_input, getArg("inputFile"));

  return;
}

template <class FieldType>
LinearParty<FieldType>::
~LinearParty(){
  delete _field;
}



#endif /* LINEARPARTY_H_ */
