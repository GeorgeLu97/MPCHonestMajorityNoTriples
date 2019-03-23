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
bool reconstruct(const vector<FieldType> code, // input
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

#endif /* BAPARTY_H_ */
