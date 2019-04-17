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
#include "ECC.h"
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/infra/Common.hpp>
// #include <libscapi/include/primitives/Prg.hpp>
// #include "HashEncrypt.h"
#include <emmintrin.h>
#include <thread>

using namespace std;
using namespace std::chrono;

// Questions:
// 1. Is ``run in parallel'' different from ``run in batch''?
// 2. count ``received bits'' includes my own bit?

// Functionality 1: concensus protocol by BGP92
// Functionality 2: broadcast protocol by CW92
// -- ^^ Currently not implemented. Assuming a broadcast channel
template <class FieldType>
class BAParty{
  
private:

  // -------- protocol properties --------
  int _myId;                    // my party id
  int _nThread;
  int _basePartyCount = 3;
  int _nActiveParties;
  vector<bool> _activeMask;     // mask to active players
  vector< shared_ptr<ProtocolPartyData> > _partySocket; // for communication
  vector<FieldType> _alpha;  // <-- the commonly known elements to all
  vector<FieldType> _beta;   //     ..one for each party
  HIM<FieldType> _peMatrix;  // common HIM matrix for PE_Broadcast
  
  // -------- private helper functions --------
  // used for debugging
  void printActiveParties();
  void printCommitee(vector<char>& commitee_mask, int myCommitee);

  // worker functions for multi-threading
  void exchangeBitWorker(bool sendBit, vector<bool>& recvBits,
                         int threadId, int nThread);
  void exchangeMsgWorker(vector<byte>& sendMsg,
                         vector< vector<byte> >& recvMsgs,
                         int msgSize, int threadId, int nThread);
  void spreadBitWorker(bool sendBit, int threadId, int nThread);
  void spreadMsgWorker(vector<byte>& msg, int threadId, int nThread);
  void gatherBitWorker(vector<bool>& recvBit, int threadId, int nThread);
  void gatherMsgWorker(vector< vector<byte> >& recvMsgs,
                       int msgSize, int threadId, int nThread);
  void scatterMsgWorker(vector< vector<byte> >& sendMsgs,
                        int msgSize, int threadId, int nThread);

    // used in base Phase King
  bool spreadBit(bool sendBit, int kingId);
  void exchangeBit(bool sendBit, vector<bool>& recvBits);
  bool universal_rounds(bool b, int* D,
                        vector<bool>& buffer, vector<bool>& buffer2);

  // used in Recursive Phase King
  int splitCommitee(vector<char>& commitee_mask, int* QCount);
  void mapCommitee(const vector<char>& commitee_mask, const int* QCount,
                    int myCommitee);
  void unmapCommitee(const vector<char>& commitee_mask, const int* QCount,
                      int myCommitee);
  void commiteeSendBit(const vector<char>& commitee_mask, const int* QCount,
                       int myCommitee, bool newb);
  void commiteeRecvBits(const vector<char>& commitee_mask, const int* QCount,
                        int myCommitee, vector<bool>& receivedBits);


public:
  // constructor and destructors
  BAParty();
  ~BAParty();
  
  // set protocol properties: (should be known before running!)
  void setParties(vector< shared_ptr<ProtocolPartyData> >& parties, int myId);
  void setHIM(HIM<FieldType>& common_HIM);
  void setAlphaBeta(vector<FieldType>& alpha, vector<FieldType>& beta);
  void setNumThreads(int numThreads);
  void getRemainingParties(vector< shared_ptr<ProtocolPartyData> >& parties);


  // -------- public functionalities --------
  // from BGP92: consensus() only on a bit (happiness).
  bool consensus_base(bool b);
  bool consensus(bool b);
  // trivial from consensus(): first dealer spread bit, then run consensus().
  void broadcastBit(bool& b, int rootId);
  void broadcastMsg(vector<byte>& msg, int msgSize, int rootId);
  void spreadMsg(vector<byte>& msg, int msgSize, int rootId);
  void scatterMsg(vector< vector<byte> >& sendMsgs, vector<byte>& recvMsg,
                  int msgSize, int rootId);
  void gatherMsg(vector<byte>& sendMsg, vector< vector<byte> >& recvMsgs,
                 int msgSize, int rootId);
  void broadcastMsgForAll(vector<byte>& sendMsg, // input
                          vector< vector<byte> >& recvMsgs,
                          int msgSize);

};


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
  
  // set default masks
  _activeMask.clear();
  _activeMask.resize(nParties, true);
  _nActiveParties = nParties +1; // counting myself
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
                  int threadId, int nThread){

  int nParties = _partySocket.size();
  char sendByte = sendBit ? 1 : 0;
  vector<char> recvBytes(nParties);

  // communicate with other parties
  for(int i=threadId; i<nParties; i+=nThread){
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
  for(int i=threadId; i<nParties; i+=nThread){
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

// send bit to all active parties
template <class FieldType>
void BAParty<FieldType>::
spreadBitWorker(bool sendBit,
                int threadId, int nThread){

  int nParties = _partySocket.size();
  char sendByte = sendBit ? 1 : 0;

  // communicate with other parties
  for(int i=threadId; i<nParties; i+=nThread){
    if(!_activeMask[i]){
      // skip inactive (eliminated) parties
      continue;
    }
    _partySocket[i]->getChannel()->write((byte*) &sendByte, 1 );
  }
  return;
}

// gather bit from all active parties
template <class FieldType>
void BAParty<FieldType>::
gatherBitWorker(vector<bool>& recvBits,
                int threadId, int nThread){

  int nParties = _partySocket.size();
  vector<char> recvBytes(nParties);
  // communicate with other parties
  for(int i=threadId; i<nParties; i+=nThread){
    if(!_activeMask[i]){
      // skip inactive (eliminated) parties
      continue;
    }

    _partySocket[i]->getChannel()->read((byte*) &(recvBytes[i]), 1 );
  }

  // convert messages into bits
  for(int i=threadId; i<nParties; i+=nThread){
    if(!_activeMask[i]){
      continue;
    }
    recvBits[i] = (recvBytes[i] == 1);
  }
  return;
}


template <class FieldType>
void BAParty<FieldType>::
scatterMsgWorker(vector< vector<byte> >& sendMsgs,
                 int msgSize, int threadId, int nThread){

  int nParties = _partySocket.size();

  // communicate with other parties
  for(int i=threadId; i<nParties; i+=nThread){
    if(!_activeMask[i]){
      // skip inactive (eliminated) parties
      continue;
    }
    sendMsgs[i].resize(msgSize);
    _partySocket[i]->getChannel()->write(sendMsgs[i].data(), msgSize);
  }

  return;
}

template <class FieldType>
bool BAParty<FieldType>::
spreadBit(bool sendBit, int kingId){
  int nParties = _partySocket.size();
  if(_myId == kingId){
    // king: send to everyone
    vector<thread> threads(_nThread);
    for(int i=0; i<_nThread; i++){
      threads[i] = thread(&BAParty::spreadBitWorker, this,
                          sendBit, i, _nThread);
    }

    for(int i=0; i<_nThread; i++){
      threads[i].join();
    }
    return sendBit;
  }else{
    // rest: ignore sendBit. Receive from king
    char recvByte;
    for(int i=0; i<nParties; i++){
      if(_activeMask[i] && _partySocket[i]->getID() == kingId){
        _partySocket[i]->getChannel()->read((byte*) &recvByte, 1);
      }
    }
    return (recvByte == 1);
  }
  
}

template <class FieldType>
void BAParty<FieldType>::
printActiveParties(){
  int nParties = _partySocket.size();
  cout << "currently active parties are: ";
  for(int i=0; i<nParties; i++){
    if(_activeMask[i]){
      cout << _partySocket[i]->getID() << " ";
    }
  }
  cout << endl;
  return;
}

template <class FieldType>
void BAParty<FieldType>::
commiteeSendBit(const vector<char>& commitee_mask, const int* QCount,
                int myCommitee, bool sendBit){
  // Send to all of the other commitee
  // -- set my commitee to in-active
  // -- broad cast to all active parties
  // -- reset my commitee to active
  mapCommitee(commitee_mask, QCount, 1-myCommitee);
  //  printActiveParties();
  vector<thread> threads(_nThread);
  for(int i=0; i<_nThread; i++){
    threads[i] = thread(&BAParty::spreadBitWorker, this,
                        sendBit, i, _nThread);
  }
  for(int i=0; i<_nThread; i++){
    threads[i].join();
  }
  unmapCommitee(commitee_mask, QCount, 1-myCommitee);
  return;
}

template <class FieldType>
void BAParty<FieldType>::
commiteeRecvBits(const vector<char>& commitee_mask, const int* QCount,
                 int myCommitee, vector<bool>& buffer){
  // Recv from all of the other commitee
  // -- set my commitee to in-active
  // -- recv from all active parties
  // -- reset my commitee to active
  mapCommitee(commitee_mask, QCount, 1-myCommitee);
  vector<thread> threads(_nThread);
  for(int i=0; i<_nThread; i++){
    threads[i] = thread(&BAParty::gatherBitWorker, this,
                        ref(buffer), i, _nThread);
  }
  for(int i=0; i<_nThread; i++){
    threads[i].join();
  }
  unmapCommitee(commitee_mask, QCount, 1-myCommitee);
  return;
}


// from BGP92: Phase King 
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
universal_rounds(bool b, int* D,
                 vector<bool>& buffer, vector<bool>& buffer2){
  int smallT = (_nActiveParties-1) /3;
  int nParties = _partySocket.size();
  int C[] = {0, 0};

  // first universal excange round.
  exchangeBit(b, buffer);
  for(int j=0; j<nParties; j++){
    if(_activeMask[j]){
      C[0] += ( !buffer[j] ); // bit == false(0)
      C[1] += ( buffer[j] );  // bit == true(1) 
    }
  }
  // second universal exchange round
  exchangeBit(C[0] >= (_nActiveParties - smallT), buffer);
  exchangeBit(C[1] >= (_nActiveParties - smallT), buffer2);
  for(int j=0; j<nParties; j++){
    if(_activeMask[j]){
      D[0] += (buffer[j]);
      D[1] += (buffer2[j]);
    }
  }
  return (D[1] > smallT);
}


template <class FieldType>
bool BAParty<FieldType>::
consensus_base(bool b){
  int smallT = (_nActiveParties-1) /3;
  int nRounds = smallT + 1;
  int nParties = _partySocket.size();
  int kingIdx = 0;
  int D[] = {0, 0};
  vector<bool> receivedBits(nParties);
  vector<bool> receivedBits2(nParties);

  bool firstGreater = true;
  
  for(int i=0; i<nRounds; i++){
    cout << "-------- base round " << i << endl;
    // skip to next active party
    while( !_activeMask[kingIdx] ){ kingIdx++; }
    int kingId;
    bool isKing;
    if(_partySocket[kingIdx]->getID() > _myId && firstGreater){
      kingId = _myId;
      isKing = true;
      firstGreater = false;
    }else{
      kingId = _partySocket[kingIdx]->getID();
      isKing = false;
    }

    // the two universal rounds
    b = universal_rounds(b, D, receivedBits, receivedBits2);
    // last king broadcast round
    bool newb = spreadBit(b, kingId);
    if(!isKing){
      int nSamePlayers = b ? D[1] : D[0];
      b = (nSamePlayers < _nActiveParties - smallT) ? newb : b;
    }// else b == newb. do nothing
  }
  
  return b;
}

template <class FieldType>
int BAParty<FieldType>::
splitCommitee(vector<char>& commitee_mask, int* QCount){
  int curCommitee = 0;
  int myCommitee; // <-- to be set in the loop
  int nParties = _partySocket.size();
  bool firstGreater = true;
  for(int i=0; i<nParties; i++){
    if(!_activeMask[i]){
      continue;
    }
    if( _myId < _partySocket[i]->getID() && firstGreater){
      // just went past myself!
      myCommitee = curCommitee;
      QCount[curCommitee]++;
      curCommitee = 1-curCommitee;
      firstGreater = false;
    }
    commitee_mask[i]= curCommitee;
    QCount[curCommitee]++;
    curCommitee = 1-curCommitee;
  }

  if(firstGreater){
    // I'm the largest ID
    myCommitee = curCommitee;
  }

  return myCommitee;
}

template <class FieldType>
void BAParty<FieldType>::
mapCommitee(const vector<char>& commitee_mask,
             const int* QCount, int myCommitee){
  int nParties = _partySocket.size();
  for(int i=0; i<nParties; i++){
    if(commitee_mask[i] == 1-myCommitee){
      _activeMask[i] = false;
    }
  }
  _nActiveParties -= QCount[1-myCommitee];
  return;
}

template <class FieldType>
void BAParty<FieldType>::
unmapCommitee(const vector<char>& commitee_mask,
               const int* QCount, int myCommitee){
  int nParties = _partySocket.size();
  for(int i=0; i<nParties; i++){
    if(commitee_mask[i] == 1-myCommitee){
      _activeMask[i] = true;
    }
  }
  _nActiveParties += QCount[1-myCommitee];  
  return;
}

template <class FieldType>
void BAParty<FieldType>::
printCommitee(vector<char>& commitee_mask, int myCommitee){
  int nParties = _partySocket.size();
  cout << "my commitee (" << myCommitee << ")  has: " << _myId;
  for(int i=0; i<nParties; i++){
    if(commitee_mask[i] == myCommitee){
      cout << " " << _partySocket[i]->getID();
    }
  }
  cout << endl;
  return;
}


// from BGP92: Recursive Phase King (RPK)
// -- with Phase King as base
// -- NOTE: assume 'parties' array is in ID order
template <class FieldType>
bool BAParty<FieldType>::
consensus(bool b){
  int sizeQw = _nActiveParties;
  // base case
  if(sizeQw <= _basePartyCount){
    return consensus_base(b);
  }

  int nParties = _partySocket.size();
  int QCount[] = {0, 0};
  int D[] = {0, 0};
  int tw = (sizeQw -1) /3;
  
  // active parties = current commitee
  // -- further divide into 2 committes, marked 1 or 0. Others marked as 2
  // -- divide commitees by consecutive active members
  vector<char> commitee_mask(nParties, 2);
  int myCommitee = splitCommitee(commitee_mask, QCount);
  // printCommitee(commitee_mask, myCommitee);

  // actual protocol
  vector<bool> receivedBits(nParties);
  vector<bool> receivedBits2(nParties);
  bool newb;
  for(int k=0; k<2; k++){
    // the two universal rounds
    cout << "-------- round " << k << endl;
    
    b = universal_rounds(b, D, receivedBits, receivedBits2);

    if(myCommitee == k){
      // if it's my commitee's turn
      // -- set the other commitee to in-active
      // -- recurse, with only my commitee being active
      // -- reset the other commitee to active
      mapCommitee(commitee_mask, QCount, myCommitee);
      newb = consensus(b);
      unmapCommitee(commitee_mask, QCount, myCommitee);

      commiteeSendBit(commitee_mask, QCount, myCommitee, newb);
    }else{
      // if not my commitee's turn
      // -- receive from the other commitee!
      // -- count most frequent value as newb
      commiteeRecvBits(commitee_mask, QCount, myCommitee, receivedBits);
      int tmpCount[] = {0, 0};
      for(int i=0; i<nParties; i++){
        if(commitee_mask[i] == 1-myCommitee){
          tmpCount[0] += (!receivedBits[i]) ? 1 : 0;
          tmpCount[1] += receivedBits[i] ? 1 : 0;
        }
      }
      newb = (tmpCount[1] > tmpCount[0]) ? 1 : 0;
    }

    if(D[b] < _nActiveParties - tw){
      b = newb;
    }
  }
  return b;
}

template <class FieldType>
void BAParty<FieldType>::
broadcastBit(bool &b, int rootId){
  b = spreadBit(b, rootId);
  b = consensus(b);
  return;
}


template <class FieldType>
void BAParty<FieldType>::
spreadMsgWorker(vector<byte>& msg,
                int threadId, int nThread) {
  int nParties = _partySocket.size();
  for (int i = threadId; i < nParties; i += nThread) {
    if (!_activeMask[i]) {
      continue;
    }
    _partySocket[i]->getChannel()->write(msg.data(), msg.size());
  }
}

template <class FieldType>
void BAParty<FieldType>::
exchangeMsgWorker(vector<byte>& sendMsg,
                  vector< vector<byte> >& recvMsgs,
                  int msgSize, int threadId, int nThread){

  int nParties = _partySocket.size();
  recvMsgs.clear();
  recvMsgs.resize(nParties, vector<byte>(msgSize));

  for (int i = threadId; i < nParties; i += nThread) {
    if (!_activeMask[i]) {
      continue;
    }

    if (_myId < _partySocket[i]->getID()) {
      // write before read
      _partySocket[i]->getChannel()->write(sendMsg.data(), msgSize);
      _partySocket[i]->getChannel()->read(recvMsgs[i].data(), msgSize);
    } else {
      // read before write
      _partySocket[i]->getChannel()->read(recvMsgs[i].data(), msgSize);
      _partySocket[i]->getChannel()->write(sendMsg.data(), msgSize);
    }
  }

  return;
}


// TODO: real implementation: broadcastForP() in BTH08 Appendix
// -- Note: will need to also implement peBroadcast()
// -- currently a stub: simply send all-to-all
template <class FieldType>
void BAParty<FieldType>::
broadcastMsgForAll(vector<byte>& sendMsg, // input
                   vector< vector<byte> >& recvMsgs, int msgSize){
  
  vector<thread> threads(_nThread);
  sendMsg.resize(msgSize);
  for(int i=0; i<_nThread; i++){
    threads[i] = thread(&BAParty::exchangeMsgWorker, this,
                        ref(sendMsg), ref(recvMsgs), msgSize, i, _nThread);
  }

  for(int i=0; i<_nThread; i++){
    threads[i].join();
  }
  return;
}


// TODO: read implementation: CW92 or check other paper
// -- currently a stub: simply send to all
template <class FieldType>
void BAParty<FieldType>::
broadcastMsg(vector<byte>& msg, int msgSize, int rootId){
  spreadMsg(msg, msgSize, rootId);
  return;
}

template <class FieldType>
void BAParty<FieldType>::
spreadMsg(vector<byte>& msg, int msgSize, int rootId){
  int nParties = _partySocket.size();

  if (rootId == _myId) {
    // I'm the sender
    msg.resize(msgSize);
    vector<thread> threads(_nThread);
    for (int i = 0; i < _nThread; i++) {
      threads[i] = thread(&BAParty::spreadMsgWorker, this,
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
          _partySocket[i]->getID() != rootId) {
        continue;
      }
      _partySocket[i]->getChannel()->read(msg.data(), msg.size());
    }
  }
  return;
}


template <class FieldType>
void BAParty<FieldType>::
scatterMsg(vector< vector<byte> >& sendMsgs,
           vector<byte>& recvMsg, int msgSize, int rootId){
  int nParties = _partySocket.size();
  
  if (rootId == _myId) {
    // I'm the sender
    sendMsgs.resize(nParties);
    vector<thread> threads(_nThread);
    for (int i = 0; i < _nThread; i++) {
      threads[i] = thread(&BAParty::scatterMsgWorker, this,
                          ref(sendMsgs), msgSize, i, _nThread);
    }
    for (int i = 0; i < _nThread; i++) {
      threads[i].join();
    }

  } else {
    // I'm a receiver
    recvMsg.resize(msgSize);
    for (int i = 0; i < nParties; i++) {
      if (!_activeMask[i] ||
          _partySocket[i]->getID() != rootId) {
        continue;
      }
      _partySocket[i]->getChannel()->read(recvMsg.data(), msgSize);
    }
  }
  return;
}


template <class FieldType>
void BAParty<FieldType>::
gatherMsgWorker(vector< vector<byte> >& recvMsgs,
                int msgSize, int threadId, int nThread){
  int nParties = _partySocket.size();

  // communicate with other parties
  for(int i=threadId; i<nParties; i+=nThread){
    if(!_activeMask[i]){
      // skip inactive (eliminated) parties
      continue;
    }
    recvMsgs[i].resize(msgSize);
    _partySocket[i]->getChannel()->read(recvMsgs[i].data(), msgSize);
  }

  return;
}


template <class FieldType>
void BAParty<FieldType>::
gatherMsg(vector<byte>& sendMsg, vector< vector<byte> >& recvMsgs,
          int msgSize, int rootId){

  int nParties = _partySocket.size();

  if(rootId == _myId){

    recvMsgs.resize(nParties);
    vector<thread> threads(_nThread);
    for(int i=0; i<_nThread; i++){
      threads[i] = thread(&BAParty::gatherMsgWorker, this,
                          ref(recvMsgs), msgSize, i, _nThread);
    }
    for(int i=0; i<_nThread; i++){
      threads[i].join();
    }
    
  }else{
    sendMsg.resize(msgSize);
    for (int i = 0; i < nParties; i++) {
      if (!_activeMask[i] ||
          _partySocket[i]->getID() != rootId) {
        continue;
      }
      _partySocket[i]->getChannel()->write(sendMsg.data(), msgSize);
    }    
  }

  return;
}

#endif /* BAPARTY_H_ */
