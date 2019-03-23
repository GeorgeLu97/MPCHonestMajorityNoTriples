#include "BAParty.h"

template <class FieldType>
bool reconstruct(const vector<FieldType> code, // input
                 vector<FieldType>& polynomial){

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
void BAParty<FieldType>::robustBatchBroadcast(vector<FieldType>& elems, bool isDealer){
  return;
}
