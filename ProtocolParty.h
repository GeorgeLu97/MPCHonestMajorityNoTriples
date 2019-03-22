#ifndef PROTOCOLPARTY_H_
#define PROTOCOLPARTY_H_

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
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/infra/Common.hpp>
#include <libscapi/include/primitives/Prg.hpp>
#include "HashEncrypt.h"
#include <emmintrin.h>
#include <thread>
#include <algorithm>
#include <utility>

#define flag_print true
#define flag_print_timings true
#define flag_print_output true


using namespace std;
using namespace std::chrono;


// Questions:
// 1. Is ``run in parallel'' different from ``run in batch''?
// 2. A consensus protocol?

template <class FieldType>
class BAParty : public Protocol, HonestMajority, Multipart{
  // from Appendex A. in BTH paper:
  // PE_Broadcast (each player spreads 1 value)
  // -- Every party P_i sends x_i to every P_j
  // -- Every party P_i computes (x''_1, ..., x''_n) = M(x'_1, ..., x'_n)
  // -- Every party P_i sends x''_i to evey P_j
  // -- Every party P_i checks if the received values are equal
  // -- Every party outputs received x_1, ..., x_n
  
  // BroadcastForP (each player spreads l values)
  // - Repeat t times (k = 0, ..., t-1):
  // -- Every party P_i sets a happy_bit = 1.
  // -- Run PE_Broadcast ceil(l/t) times, at a time
  // -- Every P_i sends to every P_j its happy_bit.
  // -- Every P_i ``AND'' all received happy bits.
  // -- Run concensus() to decide if there are unhappy players.
  // -- Run fault_localization()
  
}

template <class FieldType>
class ProtocolParty : public Protocol, public HonestMajority, MultiParty{

private:

    /**
     * N - number of parties
     * M - number of gates
     * T - number of malicious
     */

    int N, M, T, m_partyId;
    int N1, T1;
    int bigT;
    int times; //number of times to run the run function
    int iteration; //number of the current iteration

    Measurement* timer;
    TemplateField<FieldType> *field;
    vector<shared_ptr<ProtocolPartyData>> parties;

    vector<shared_ptr<ProtocolPartyData>> activeParties; /* partyId -> comm ptr */
    vector<int> activePartyIDs; //size N1
    bool selfActive;

    vector<tuple<FieldType,FieldType,FieldType>> randomTriplesArray;

    vector<byte> h;//a string accumulated that should be hashed in the comparing views function.

    ProtocolTimer* protocolTimer;
    int currentCirciutLayer = 0;
    int randomSharesOffset = 0;

    string s;
    int numOfInputGates, numOfOutputGates;
    string inputsFile, outputFile;

    /* TODO: probably get rid of all of these for now, 
    only useful for optimization with fixed parties */


    HIM<FieldType> m;

    boost::asio::io_service io_service;
    ArithmeticCircuit circuit;
    vector<FieldType> gateValueArr; // the value of the gate (for my input and output gates)
    vector<FieldType> gateShareArr; // my share of the gate (for all gates)
    //vector<FieldType> alpha; // N distinct non-zero field elements

    vector<long> myInputs;

    void trimZeroes(vector<FieldType>& a);
    void polyAdd2(vector<FieldType>& a, vector<FieldType>& b, vector<FieldType>& sum);
    void polyMult2(vector<FieldType>& a, vector<FieldType>& b, vector<FieldType>& product);
    void polyAdd(vector<vector<FieldType>>& addends, vector<FieldType>& sum);
    void polyMult(vector<vector<FieldType>>& multiplicands, vector<FieldType>& product);
    void scalMult(FieldType scalar, vector<FieldType>& poly, vector<FieldType>& product);
    vector<byte> encodePair1(pair<int,int> f);
    pair<int,int> decodePair1(vector<byte>& f);

public:
//    ProtocolParty(int n, int id,string fieldType, string inputsFile, string outputFile, string circuitFile,
//             int groupID = 0);
    ProtocolParty(int argc, char* argv[]);

    //lagrange interpolation
    //TODO: make robust?
    void lagrange(vector<FieldType>& alpha, vector<FieldType>& x, vector<FieldType>& coeffs);

    void roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round);
    void exchangeData(vector<vector<byte>> &sendBufs,vector<vector<byte>> &recBufs, int first, int last);
    void roundFunctionSyncBroadcast(vector<byte> &message, vector<vector<byte>> &recBufs);

    void roundFunctionSyncForP1(vector<byte>& selfMessage, vector<vector<byte>> &recBufs);
    void recDataToP1(vector<byte>& selfMessage, vector<vector<byte>> &recBufs, int first, int last);

    void sendDataFromP1(vector<byte> &sendBuf, int first, int last);
    void sendFromP1(vector<byte> &sendBuf);

    void decodeFieldElts(vector<byte>& input, vector<FieldType>& output);
    void encodeFieldElts(vector<FieldType>& input, vector<byte>& output);

    void printVector(vector<FieldType>& input);

    int counter = 0;

    /**
     * This method runs the protocol:
     * 1. Preparation Phase
     * 2. Input Phase
     * 3. Computation Phase
     * 4. Verification Phase
     * 5. Output Phase
     */
    void run() override;

    bool hasOffline() {
        return true;
    }


    bool hasOnline() override {
        return true;
    }

    /**
     * This method runs the protocol:
     * Preparation Phase
     */
    void runOffline() override;

    /**
     * This method runs the protocol:
     * Input Phase
     * Computation Phase
     * Verification Phase
     * Output Phase
     */
    void runOnline() override;

    /**
     * This method reads text file and inits a vector of Inputs according to the file.
     */
    void readMyInputs();

    /**
     * We describe the protocol initialization.
     * In particular, some global variables are declared and initialized.
     */
    void initializationPhase();


    /**
     * A random double-sharing is a pair of two sharings of the same random value, where the one sharing is
     * of degree t, and the other sharing is of degree 2t. Such random double-sharing are of big help in the
     * multiplication protocol.
     * We use hyper-invertible matrices to generate random double-sharings. The basic idea is as follows:
     * Every party generates one random double-sharing. These n double-sharings are processes through a
     * hyper-invertible matrix. From the resulting n double-sharings, t are checked to be valid (correct degree,
     * same secret), and t are then kept as “good” double-sharings. This is secure due to the diversion property
     * of hyper-invertible matrix: We know that n − t of the input double-sharings are good. So, if there are t
     * valid output double-sharings, then all double-sharings must be valid. Furthermore, the adversary knows
     * his own up to t input double-sharings, and learns t output double sharings. So, n − 2t output double
     * sharings are random and unknown to the adversary.
     * For the sake of efficiency, we do not publicly reconstruct t of the output double-sharings. Rather, we
     * reconstruct 2t output double sharings, each to one dedicated party only. At least t of these parties are
     * honest and correctly validate the reconstructed double-sharing.
     *
     * The goal of this phase is to generate “enough” double-sharings to evaluate the circuit. The double-
     * sharings are stored in a buffer SharingBuf , where alternating a degree-t and a degree-2t sharing (of the same secret)
     * is stored (more precisely, a share of each such corresponding sharings is stored).
     * The creation of double-sharings is:
     *
     * Protocol Generate-Double-Sharings:
     * 1. ∀i: Pi selects random value x-(i) and computes degree-t shares x1-(i) and degree-2t shares x2-(i).
     * 2. ∀i,j: Pi sends the shares x1,j and X2,j to party Pj.
     * 3. ∀j: Pj applies a hyper-invertible matrix M on the received shares, i.e:
     *      (y1,j,..., y1,j) = M(x1,j,...,x1,j)
     *      (y2,j,...,y2,j) = M (x2,j,...,x2,)
     * 4. ∀j, ∀k ≤ 2t: Pj sends y1,j and y2,j to Pk.
     * 5. ∀k ≤ 2t: Pk checks:
     *      • that the received shares (y1,1,...,y1,n) are t-consistent,
     *      • that the received shares (y2,1,...,y2,n) are 2t-consistent, and
     *      • that both sharings interpolate to the same secret.
     *
     * We use this algorithm, but extend it to capture an arbitrary number of double-sharings.
     * This is, as usual, achieved by processing multiple buckets in parallel.
     */
    bool preparationPhase();


    /**
     * The input phase proceeds in two steps:
     * First, for each input gate, the party owning the input creates shares for that input by choosing a random coefficients for the polynomial
     * Then, all the shares are sent to the relevant party
     */
    void inputPhase();

    void generateRandomShares(int numOfRandoms, vector<FieldType> &randomElementsToFill);
    void getRandomShares(int numOfRandoms, vector<FieldType> &randomElementsToFill);
    void generateRandomSharesWithCheck(int numOfRnadoms, vector<FieldType>& randomElementsToFill);

    bool doubleShareRandom(int degree1, int degree2, vector<tuple<FieldType, FieldType>>& randomDoubles);
    bool doubleShareRandomVerifyOne(int degree1, int degree2, 
    vector<FieldType>& d1coefficients, vector<FieldType>& d2coefficients, int partyID, 
    vector<vector<vector<vector<byte>>>>& inputs, vector<vector<vector<vector<byte>>>>& outputs);
    pair<int,int> doubleShareRandomVerifyAll(int degree1, int degree2, 
    vector<FieldType>& d1coefficients, vector<FieldType>& d2coefficients,
    vector<vector<vector<byte>>>& inputs, vector<vector<vector<byte>>>& outputs);

    //modification - do fault detection/localization here.
    void generateTriples(vector<tuple<FieldType, FieldType, FieldType>>& randomTriples);

    void removeParties(int party1, int party2);

    pair<int,int> IOConsistency(vector<vector<vector<byte>>>& sentIO, vector<vector<vector<byte>>>& receivedIO);
    pair<int,int> JudgeIOConsistency(vector<vector<byte>>& mySend, vector<vector<byte>>& myReceive, 
    vector<vector<vector<byte>>>& allSends, vector<vector<vector<byte>>>& allRecs);
    pair<int,int> NonJudgeIOConsistency(vector<vector<byte>>& mySend, vector<vector<byte>>& myReceive);



    //return happiness
    //uses P'
    //takes n-2t shares and publicly reconstructs or fails w/ unhappiness
    bool reconstructPublic(vector<FieldType>& secretShares, int degree, vector<FieldType>& reconstructedSecrets);
    bool reconstructPublicVerifyOne(vector<FieldType>& secretShares, int degree, int partyID, 
        vector<vector<vector<vector<byte>>>>& inputs, vector<vector<vector<vector<byte>>>>& outputs);

    pair<int,int> reconstructPublicVerifyAll(vector<FieldType>& secretShares, int degree, 
        vector<vector<vector<byte>>>& inputs, vector<vector<vector<byte>>>& outputs);



    //trusted party
    //partyID can be non active party
    FieldType reconstructPrivate(FieldType secretShare, int degree, int partyID);

    /**
     * Check whether given points lie on polynomial of degree d.
     * This check is performed by interpolating x on the first d + 1 positions of α and check the remaining positions.
     */
    bool checkConsistency(vector<FieldType>& x, int d);

    bool checkConsistency(vector<FieldType>& p, vector<FieldType>& x, int d);

    void openShare(int numOfRandomShares, vector<FieldType> &Shares, vector<FieldType> &secrets);

    FieldType evaluatePolynomial(FieldType x, vector<FieldType>& coefficients, int start, int degree);

    FieldType evaluatePolynomial(FieldType x, vector<FieldType>& coefficients);

    bool isHappy(bool selfHappiness);

    /**
     * Process all multiplications which are ready.
     * Return number of processed gates.
     */
    int processMultiplications(int lastMultGate);

    int processNotMult();

    /**
     * Walk through the circuit and evaluate the gates. Always take as many gates at once as possible,
     * i.e., all gates whose inputs are ready.
     * We first process all random gates, then alternately process addition and multiplication gates.
     */
    void computationPhase();

    /**
     * The cheap way: Create a HIM from the αi’s onto ZERO (this is actually a row vector), and multiply
     * this HIM with the given x-vector (this is actually a scalar product).
     * The first (and only) element of the output vector is the secret.
     */
    FieldType interpolate(vector<FieldType>& x, int d);

    bool interpolateCoefficients(vector<FieldType>& x, int d, vector<FieldType>& coeffs);

    vector<FieldType> getAlpha();

    void interpolateFull(vector<FieldType>& alpha, vector<FieldType>& x, 
        vector<FieldType>& beta, vector<FieldType>& result);

    bool interpolateCoefficientsFull(vector<FieldType>& alpha, vector<FieldType>& x, vector<FieldType>& coeffs);

    /**
     * Walk through the circuit and reconstruct output gates.
     */
    void outputPhase();

    ~ProtocolParty();



};

template <class FieldType> 
void ProtocolParty<FieldType>::printVector(vector<FieldType>& input) {
    if(input.size() == 0) {return;}
    cout << "[" << input[0];
    for(int i = 1; i < input.size(); i++) {
        cout << ", " << input[i];
    }
    cout << "]" << endl;
}

template <class FieldType>
ProtocolParty<FieldType>::ProtocolParty(int argc, char* argv[]) : Protocol("PerfectlySecureLinearCommunication", argc, argv)
{
    cout << "Start Constructor" << endl;
    
    string circuitFile = this->getParser().getValueByKey(arguments, "circuitFile");
    this->times = stoi(this->getParser().getValueByKey(arguments, "internalIterationsNumber"));
    string fieldType = this->getParser().getValueByKey(arguments, "fieldType");
    m_partyId = stoi(this->getParser().getValueByKey(arguments, "partyID"));
    int n = stoi(this->getParser().getValueByKey(arguments, "partiesNumber"));
    string outputTimerFileName = circuitFile + "Times" + to_string(m_partyId) + fieldType + ".csv";
    ProtocolTimer p(times, outputTimerFileName);

    this->protocolTimer = new ProtocolTimer(times, outputTimerFileName);

    vector<string> subTaskNames{"Offline", "preparationPhase", "Online", "inputPhase", "ComputePhase", "outputPhase"};
    timer = new Measurement(*this, subTaskNames);

    cout << "Fin Parsing" << endl;

    if(fieldType.compare("ZpMersenne") == 0) {
        field = new TemplateField<FieldType>(2147483647);
    } else if(fieldType.compare("ZpMersenne61") == 0) {
        field = new TemplateField<FieldType>(0);
    } else if(fieldType.compare("ZpKaratsuba") == 0) {
        field = new TemplateField<FieldType>(0);
    } else if(fieldType.compare("GF2E") == 0) {
        field = new TemplateField<FieldType>(8);
    } else if(fieldType.compare("Zp") == 0) {
        field = new TemplateField<FieldType>(2147483647);
    }


    N = n; //at least 3
    T = n/3; //at least 1 
    N1 = N;
    T1 = T;
    bigT = N - 2 * T;
    this->inputsFile = this->getParser().getValueByKey(arguments, "inputFile");
    this->outputFile = this->getParser().getValueByKey(arguments, "outputFile");

    s = to_string(m_partyId);
    circuit.readCircuit(circuitFile.c_str());
    circuit.reArrangeCircuit();
    M = circuit.getNrOfGates();
    numOfInputGates = circuit.getNrOfInputGates();
    numOfOutputGates = circuit.getNrOfOutputGates();
    myInputs.resize(numOfInputGates);
    counter = 0;


    //comm->ConnectionToServer(s);

    //boost::asio::io_service io_service;


    cout << "TryComm" << endl;

    MPCCommunication comm;
    string partiesFile = this->getParser().getValueByKey(arguments, "partiesFile");

    parties = comm.setCommunication(io_service, m_partyId, N, partiesFile);

    selfActive = true;
    for(int i = 0; i < n; i++) {
        activePartyIDs.push_back(i);
        if(i < m_partyId) { activeParties.push_back(parties[i]); }
        else if(i == m_partyId ) { activeParties.push_back(nullptr); }
        else { activeParties.push_back(parties[i-1]); }
    }

    cout << "ActiveParties" << endl;

    string tmp = "init times";
    //cout<<"before sending any data"<<endl;
    byte tmpBytes[20];
    for (int i=0; i<parties.size(); i++){
        if (parties[i]->getID() < m_partyId){
            parties[i]->getChannel()->write(tmp);
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
        } else {
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
            parties[i]->getChannel()->write(tmp);
        }
    }

    readMyInputs();
    cout << "InputsRead" << endl;

    auto t1 = high_resolution_clock::now();
    initializationPhase(/*matrix_him, matrix_vand, m*/);

    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds initializationPhase: " << duration << endl;
    }
}


template<class FieldType> 
void ProtocolParty<FieldType>::trimZeroes(vector<FieldType>& a) {
    int i = a.size();
    FieldType Zero = *(field->GetZero());
    while(i > 0 && a[i-1] == Zero) { i--; }
    a.resize(i);
}

template <class FieldType>
void ProtocolParty<FieldType>::polyAdd2(vector<FieldType>& a, vector<FieldType>& b, vector<FieldType>& sum) {
    vector<FieldType> tempSum(a);
    if(b.size() > a.size()) { tempSum.resize(b.size()); }
    for(int i = 0; i < b.size(); i++) {
        tempSum[i] += b[i];
    }
    sum = tempSum;
    trimZeroes(sum);
}

template <class FieldType>
void ProtocolParty<FieldType>::polyMult2(vector<FieldType>& a, vector<FieldType>& b, vector<FieldType>& product) {
    vector<FieldType> tempProduct(a.size() + b.size() - 1); 
    for(int i = 0; i < a.size(); i++) {
        for(int j = 0; j < b.size(); j++) {
            tempProduct[i + j] += a[i] * b[j];
        }
    }
    product = tempProduct;
    trimZeroes(product);
}

template <class FieldType>
void ProtocolParty<FieldType>::polyAdd(vector<vector<FieldType>>& addends, vector<FieldType>& sum) {
    vector<FieldType> partialSum;
    for(auto it = addends.begin(); it != addends.end(); it++) {
        polyAdd2(partialSum, *it, partialSum);
    }
    sum = partialSum;
}

template <class FieldType>
void ProtocolParty<FieldType>::polyMult(vector<vector<FieldType>>& multiplicands, vector<FieldType>& product) {
    vector<FieldType> partialProduct;
    partialProduct.push_back(*(field->GetOne()));
    for(auto it = multiplicands.begin(); it != multiplicands.end(); it++) {
        polyMult2(partialProduct, *it, partialProduct);
    }
    product = partialProduct;
}

template <class FieldType>
void ProtocolParty<FieldType>::scalMult(FieldType scalar, vector<FieldType>& poly, vector<FieldType>& product) {
    product = poly;
    for(int i = 0; i < product.size(); i++) {
        product[i] *= scalar;
    }
}

//lagrange interpolation - from (alpha,x) to coefficients of polynomial
//alpha is x coefficient, x is "y" coefficients 
template <class FieldType>
void ProtocolParty<FieldType>::lagrange(vector<FieldType>& alpha, vector<FieldType>& x, vector<FieldType>& coeffs) {
    vector<FieldType> unit;
    unit.push_back(*(field->GetOne()));
    vector<vector<FieldType>> linearPolys(alpha.size());
    vector<vector<FieldType>> prefixPolys(alpha.size()); //prefix[i] = (x-alpha_{0 ... i-1})
    vector<vector<FieldType>> suffixPolys(alpha.size()); //suffix[i] = (x-alpha_{i + 1 ... n-1})
    vector<vector<FieldType>> exclusionPolys(alpha.size()); //exP[i] should be 0 all alpha except alpha[i]
    vector<vector<FieldType>> xExclusionPolys(alpha.size());
    prefixPolys[0] = unit;
    suffixPolys[alpha.size()-1] = unit;

    for(int i = 0; i < alpha.size(); i++) {
        linearPolys[i].push_back(*(field->GetZero()) - alpha[i]); //-alpha_i
        linearPolys[i].push_back(*(field->GetOne()));
    }
    
    for(int i = 1; i < alpha.size(); i++) {
        polyMult2(prefixPolys[i-1], linearPolys[i-1], prefixPolys[i]);
        polyMult2(suffixPolys[alpha.size() - i], linearPolys[alpha.size() - i], suffixPolys[alpha.size() - i - 1]);
    }


    /*
    for(int i = 0; i < alpha.size(); i++) {
        for(int j = 0; j < prefixPolys[i].size(); j++) {
            cout << prefixPolys[i][j] << "x^" << j << " + ";
        }
        cout << endl;
    }

    */

    for(int i = 0; i < alpha.size(); i++) {
        polyMult2(prefixPolys[i], suffixPolys[i], exclusionPolys[i]);
        FieldType s = (*(field->GetOne())) / evaluatePolynomial(alpha[i], exclusionPolys[i]);
        scalMult(s * x[i], exclusionPolys[i], xExclusionPolys[i]);
        //xExclusionPolys[i](alpha_i) = x, xExclusionPolys[i](alpha_j) = 0
    }

    polyAdd(xExclusionPolys, coeffs);

    /*
    for(int i = 0; i < alpha.size(); i++) {
        for(int j = 0; j < exclusionPolys.size(); j++) {
            cout << "xExclusionTest: " << i << " " << j << " " << alpha[i] << " " << x[i] << " " << evaluatePolynomial(alpha[i], xExclusionPolys[j]) << endl;
        }
    }

    for(int i = 0; i < alpha.size(); i++) {
        cout << "LagrangeTest: " << alpha[i] << " " << x[i] << " " << evaluatePolynomial(alpha[i], coeffs) << endl;
    } */
}

template <class FieldType>
void transpose(vector<vector<FieldType>>& input, vector<vector<FieldType>>& output) {
    if(input.size()) {
        output.resize(input[0].size());
        for(auto i = 0; i < input[0].size(); i++) {
            output[i].resize(input.size());
            for(int j = 0; j < input.size(); j++) {
                output[i][j] = input[j][i];
            }
        } 
    } else {
        output.clear();
        return;
    }
}

//sendbufs/recbufs correspond to activeParties
template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round) {

    //cout<<"in roundFunctionSync "<< round<< endl;

    int numThreads = 10;//parties.size();
    int numPartiesForEachThread;

    if (activeParties.size() <= numThreads){
        numThreads = activeParties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (activeParties.size() + numThreads - 1)/ numThreads;
    }
    
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= activeParties.size()) {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, activeParties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::exchangeData(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int first, int last){


    //cout<<"in exchangeData";
    for (int i=first; i < last; i++) {
        if(m_partyId == activePartyIDs[i]) {
            recBufs[i] = move(sendBufs[i]);
        } else if ((m_partyId) < activePartyIDs[i]) {
            //send shares to my input bits
            activeParties[i]->getChannel()->write(sendBufs[i].data(), sendBufs[i].size());
            //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;

            //receive shares from the other party and set them in the shares array
            activeParties[i]->getChannel()->read(recBufs[i].data(), recBufs[i].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
        } else{
            //receive shares from the other party and set them in the shares array
            activeParties[i]->getChannel()->read(recBufs[i].data(), recBufs[i].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;

            //send shares to my input bits
            activeParties[i]->getChannel()->write(sendBufs[i].data(), sendBufs[i].size());
            //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;
        }

    }


}

//broadcast data
template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSyncBroadcast(vector<byte> &message, vector<vector<byte>> &recBufs) {
    vector<vector<byte>> sendBufs;
    for(auto it = recBufs.begin(); it != recBufs.end(); it++) { sendBufs.push_back(message); }
    roundFunctionSync(sendBufs, recBufs, 19);
}


template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSyncForP1(vector<byte>& selfMessage, vector<vector<byte>> &recBufs) {

    //cout<<"in roundFunctionSyncBroadcast "<< endl;

    int numThreads = 10; //parties.size();
    int numPartiesForEachThread;

    if (activeParties.size() <= numThreads){
        numThreads = activeParties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (activeParties.size() + numThreads - 1)/ numThreads;
    }

    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= activeParties.size()) {
            threads[t] = thread(&ProtocolParty::recDataToP1, this,  ref(selfMessage), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::recDataToP1, this, ref(selfMessage), ref(recBufs), t * numPartiesForEachThread, activeParties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::recDataToP1(vector<byte>& selfMessage, vector<vector<byte>> &recBufs, int first, int last){
    for (int i=first; i < last; i++) {
        if(activePartyIDs[i] == m_partyId) {
            recBufs[i] = move(selfMessage);
        } else {
            activeParties[i]->getChannel()->read(recBufs[i].data(), recBufs[i].size());
        }
    }
}



template <class FieldType>
void ProtocolParty<FieldType>::sendFromP1(vector<byte> &sendBuf) {

    //cout<<"in roundFunctionSyncBroadcast "<< endl;

    int numThreads = 10; //activeParties.size();
    int numPartiesForEachThread;

    if (activeParties.size() <= numThreads){
        numThreads = activeParties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (activeParties.size() + numThreads - 1)/ numThreads;
    }

    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= activeParties.size()) {
            threads[t] = thread(&ProtocolParty::sendDataFromP1, this,  ref(sendBuf),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::sendDataFromP1, this, ref(sendBuf), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}

template <class FieldType>
void ProtocolParty<FieldType>::sendDataFromP1(vector<byte> &sendBuf, int first, int last){
    for(int i=first; i < last; i++) {
        if(activePartyIDs[i] != m_partyId) {
            activeParties[i]->getChannel()->write(sendBuf.data(), sendBuf.size());
        }
    }
}

//
template <class FieldType>
void ProtocolParty<FieldType>::encodeFieldElts(vector<FieldType>& input, vector<byte>& output){
    int fieldByteSize = field->getElementSizeInBytes();
    output.resize(fieldByteSize * input.size());
    for(int j = 0; j < input.size(); j++) {
        field->elementToBytes(output.data() + (j * fieldByteSize), input[j]);
    }
}

template <class FieldType>
void ProtocolParty<FieldType>::decodeFieldElts(vector<byte>& input, vector<FieldType>& output){
    int fieldByteSize = field->getElementSizeInBytes();
    if(input.size() % fieldByteSize) { cout << "input byte misalignment" << endl; }
    output.resize(input.size() / fieldByteSize);

    for(int j = 0; j < output.size(); j++) {
        output[j] = field->bytesToElement(input.data() + fieldByteSize * j);
    }
}

template <class FieldType>
bool ProtocolParty<FieldType>::checkConsistency(vector<FieldType>& x, int d)
{
    vector<FieldType> alpha = getAlpha();
    return checkConsistency(alpha, x, d);
}

/* p is input fieldelt's */
//TODO: Only need T + t' = n' - t' consistent values
//TODO: make this robust rather than just detectable
//Or, incorporate dispute control here?
template <class FieldType>
bool ProtocolParty<FieldType>::checkConsistency(vector<FieldType>& p, vector<FieldType>& x, int d)
{
    if(p.size() <= d+1) {
        return true;
    } else {
        vector<FieldType> alpha_until_d(d + 1);
        vector<FieldType> alpha_from_d(p.size() - 1 - d);
        vector<FieldType> x_until_d(d + 1);
        vector<FieldType> y(p.size() - 1 - d); // the result of multiplication

        for (int i = 0; i < d + 1; i++) {
            alpha_until_d[i] = p[i];
            x_until_d[i] = x[i];
        }
        for (int i = d + 1; i < N; i++) {
            alpha_from_d[i - (d + 1)] = p[i];
        }
        // Interpolate first d+1 positions of (alpha,x)
        HIM<FieldType> matrix(p.size() - 1 - d, d + 1, field); // slices, only positions from 0..d
        matrix.InitHIMByVectors(alpha_until_d, alpha_from_d);
        matrix.MatrixMult(x_until_d, y);

        // compare that the result is equal to the according positions in x
        for (int i = 0; i < y.size(); i++)   // n-d-2 or n-d-1 ??
        {
            if (y[i] != x[d + 1 + i]) {
                return false;
            }
        }
        return true;
    }
    return true;
}

// Interpolate polynomial at position Zero
//degree shouldn't matter as long as we have checked consistency.
//x must be n' long
//TODO: make robust
template <class FieldType>
FieldType ProtocolParty<FieldType>::interpolate(vector<FieldType>& x, int d)
{
    (void)d; //unused for now, but will be needed for error correction
    vector<FieldType> result(1);
    vector<FieldType> beta(1);
    vector<FieldType> alpha = getAlpha();
    beta[0] = *(field->GetZero()); 
    //takes advantage of fact that constant term of polynomial = p(0)
    interpolateFull(alpha, x, beta, result);
    return result[0];
}

//alpha should be N' long
template <class FieldType>
vector<FieldType> ProtocolParty<FieldType>::getAlpha()
{
    vector<FieldType> alpha(activePartyIDs.size());
    for(int i = 0; i < activePartyIDs.size(); i++) {
        alpha[i] = field->GetElement(activePartyIDs[i] + 1);
    }
    return alpha;
}

// Interpolate polynomial at position Zero
template <class FieldType>
void ProtocolParty<FieldType>::interpolateFull(vector<FieldType>& alpha, vector<FieldType>& x, vector<FieldType>& beta,
 vector<FieldType>& result)
{
    //vector<FieldType> y(N); // result of interpolate
    HIM<FieldType> matrix(beta.size(),alpha.size(),field);
    matrix.InitHIMByVectors(alpha, beta);
    matrix.MatrixMult(x, result);
}


//gives sum coeff[i]x^i
//where coeff = coefficients[start, start + degree]
//degree = length - 1 = last index used
template <class FieldType>
FieldType ProtocolParty<FieldType>::evaluatePolynomial(FieldType x, vector<FieldType>& coefficients, int start, int degree)
{
    FieldType res = *(field->GetZero());
    //if degree is too large, 
    int lastindex = start + degree;
    if(lastindex >= coefficients.size()) { lastindex = coefficients.size() - 1; }
    for(int i = start + degree; i >= start; i--) {
        res *= x;
        res += coefficients[i];
    }
    return res;
}

template <class FieldType>
FieldType ProtocolParty<FieldType>::evaluatePolynomial(FieldType x, vector<FieldType>& coefficients)
{
    return evaluatePolynomial(x, coefficients, 0, coefficients.size() - 1);
}

template <class FieldType> 
bool ProtocolParty<FieldType>::isHappy(bool selfHappiness) {
    vector<byte> happiness;
    vector<vector<byte>> otherHappiness(activeParties.size());
    happiness.push_back(selfHappiness ? 1 : 0);
    for(int i = 0; i < activeParties.size(); i++) {
        otherHappiness[i].resize(1);
    }
    //TODO: byzantine
    roundFunctionSyncBroadcast(happiness, otherHappiness);
    if(!selfHappiness) { return false; }
    for(int i = 0; i < activeParties.size(); i++) {
        if(otherHappiness[i][0] != 1) { cout << "UnHappy!!!!" << endl; return false; }
    }
    return true;
}


template <class FieldType>
vector<byte> ProtocolParty<FieldType>::encodePair1(pair<int,int> f) {
    vector<byte> results(8);

    unsigned int f1 = (unsigned int)f.first;
    unsigned int f2 = (unsigned int)f.second;

    //encode 
    for(int i = 0; i < sizeof(int); i++) {
        results[2*i] = (unsigned char)f1;
        results[2*i + 1] = (unsigned char)f2;
        f1 >>= 8;
        f2 >>= 8;
    }

    return results;
}

template <class FieldType>
pair<int,int> ProtocolParty<FieldType>::decodePair1(vector<byte>& f) {
    pair<int,int> g;

    for(int i = sizeof(int) - 1; i >= 0; i--) {
        g.first <<= 8;
        g.second <<= 8;
        g.first += (int)((unsigned int)f[2*i]);
        g.second += (int)((unsigned int)f[2*i + 1]);
    }
    return g;
}

//returns field elt if reconstructed party, else zero
template <class FieldType> 
FieldType ProtocolParty<FieldType>::reconstructPrivate(FieldType secretShare, int degree, int partyID) {
    vector<byte> ssb(field->getElementSizeInBytes());
    field->elementToBytes(ssb.data(), secretShare);
    if(m_partyId == partyID) {
        cout << m_partyId << "rp1" << endl;
        vector<vector<byte>> allShares(N1);
        vector<FieldType> allSharesTrue(N1);
        for(int i = 0; i < N1; i++) {
            allShares[i].resize(field->getElementSizeInBytes());
        } 
        cout << m_partyId << "rp2" << endl;
        roundFunctionSyncForP1(ssb, allShares);
        cout << m_partyId << "rp3" << endl;
        for(int i = 0; i < N1; i++) {
            allSharesTrue[i] = field->bytesToElement(allShares[i].data());
        }
        FieldType ans = interpolate(allSharesTrue, degree);

        return ans;
    } else if(selfActive) {
        cout << m_partyId << "writeSize" << ssb.size() << endl;
        activeParties[partyID]->getChannel()->write(ssb.data(), ssb.size());
        cout << m_partyId << "rp5" << endl;
        //parties[partyID]->getChannel()->read(ssb.data(), ssb.size());
        //return field->bytesToElement(ssb.data());
        return *(field->GetZero());
    }
}

//TODO: add fault localization
//we can group shares for better round efficiency or separate for less dis impact
//reconstructShare is appended to, not overwritten
template <class FieldType>
bool ProtocolParty<FieldType>::reconstructPublic(vector<FieldType>& secretShares, int degree, vector<FieldType>& reconstructedSecrets)
{

    bool happiness = true;
    vector<vector<FieldType>> ujs;
    for(auto it = activePartyIDs.begin(); it != activePartyIDs.end(); it++) {
        vector<FieldType> partyPoly;
        for(int start = 0; start < secretShares.size(); start += bigT) {
            partyPoly.push_back(evaluatePolynomial(field->GetElement((*it) + 1), secretShares, start, bigT - 1));
        }
        ujs.push_back(partyPoly);
    }

    vector<vector<byte>> sendBufs(N1);
    vector<vector<byte>> recBufs(N1);

    vector<byte> uiSend;
    vector<vector<byte>> uiRec(N1);
    vector<vector<byte>> uiSendAll(N1);

    for(int i = 0; i < activeParties.size(); i++) {
        encodeFieldElts(ujs[i], sendBufs[i]);
        recBufs[i].resize(sendBufs[i].size());
    }

    //sendBufs should have degree bigT poylnomial
    roundFunctionSync(sendBufs, recBufs, 13);
    vector<vector<FieldType>> ujd(N1);

    for(int i = 0; i < activeParties.size(); i++) {
        decodeFieldElts(recBufs[i], ujd[i]);
    }
    
    vector<vector<FieldType>> ujd_transpose;
    transpose(ujd, ujd_transpose);

    for(auto it = ujd_transpose.begin(); it != ujd_transpose.end(); it++) {
        if(!checkConsistency(*it, degree)) {
            happiness = false;
            break;
        }
    }

    if(!isHappy(happiness)) { 
        vector<vector<vector<byte>>> inputs;
        vector<vector<vector<byte>>> outputs;
        inputs.push_back(sendBufs);
        inputs.push_back(uiSendAll);
        outputs.push_back(recBufs);
        outputs.push_back(uiRec);
        reconstructPublicVerifyAll(secretShares, degree, inputs, outputs);
        return false;
    }

    //recovered ui;
    vector<FieldType> ui;
    for(auto it = ujd_transpose.begin(); it != ujd_transpose.end(); it++) {
        ui.push_back(interpolate(*it, degree));
    }

    //broadcast computed ui's
    encodeFieldElts(ui, uiSend);
    for(auto it = uiRec.begin(); it != uiRec.end(); it++) {
        it->resize(uiSend.size());
    }
    //uisend should be constant
    // = s1 + s2b + s3b^2 + ... where b = group elt and s1 is true secret. 
    roundFunctionSyncBroadcast(uiSend, uiRec);

    vector<vector<FieldType>> uj(N1);

    for(int i = 0; i < activePartyIDs.size(); i++) {
        decodeFieldElts(uiRec[i], uj[i]);
    }

    vector<vector<FieldType>> uj_transpose;
    transpose(uj, uj_transpose);
    
    //uiSend should be repeated field elt
    //uiRec should be 
    for(auto it = uj_transpose.begin(); it != uj_transpose.end(); it++) {
        if(!checkConsistency(*it, bigT - 1)) {
            happiness = false;
            break;
        }
    }

    if(!isHappy(happiness)) {
        vector<vector<vector<byte>>> inputs;
        vector<vector<vector<byte>>> outputs;
        for(int i = 0; i < N1; i++) {
            uiSendAll.push_back(uiSend);
        }
        inputs.push_back(sendBufs);
        inputs.push_back(uiSendAll);
        outputs.push_back(recBufs);
        outputs.push_back(uiRec);
        reconstructPublicVerifyAll(secretShares, degree, inputs, outputs); 
        return false; 
    }

    //ui -> si values to coefficients
    for(auto it = uj_transpose.begin(); it != uj_transpose.end(); it++) {
        vector<FieldType> alpha = getAlpha();
        vector<FieldType> si;
        lagrange(alpha, *it, si);

        si.resize(bigT);

        for(auto it2 = si.begin(); it2 != si.end(); it2++) { reconstructedSecrets.push_back(*it2); }
    }

    reconstructedSecrets.resize(secretShares.size());
    return true;
}

template <class FieldType>
bool ProtocolParty<FieldType>::reconstructPublicVerifyOne(vector<FieldType>& secretShares, int degree, int partyID,
        vector<vector<vector<vector<byte>>>>& inputs, vector<vector<vector<vector<byte>>>>& outputs)
{
    //It's possible someone lies about their secretShares, so we should check that too in verifyAll
    //inputs[0][partyId] = sendBufs
    //inputs[1][partyId] = uiSend x n

    //outputs[0][partyId] = recBufs
    //outputs[1][partyId] = uiRec
    bool happiness = true;
    vector<vector<FieldType>> ujs;
    for(auto it = activePartyIDs.begin(); it != activePartyIDs.end(); it++) {
        vector<FieldType> partyPoly;
        for(int start = 0; start < secretShares.size(); start += bigT) {
            partyPoly.push_back(evaluatePolynomial(field->GetElement((*it) + 1), secretShares, start, bigT - 1));
        }
        ujs.push_back(partyPoly);
    }

    vector<vector<byte>> sendBufs(N1);
    vector<vector<byte>> recBufs(N1);
    for(int i = 0; i < activeParties.size(); i++) {
        encodeFieldElts(ujs[i], sendBufs[i]);
        recBufs[i].resize(ujs[i].size());
    }


    //check == inputs;
    //sendBufs should have degree bigT poylnomial
    //roundFunctionSync(sendBufs, recBufs, 13);
    //first input
    if(sendBufs != inputs[0][partyID]) {
        return false;
    } else {
        recBufs = outputs[0][partyID];
    }


    vector<vector<FieldType>> ujd(N1);

    for(int i = 0; i < activeParties.size(); i++) {
        decodeFieldElts(recBufs[i], ujd[i]);
    }
    
    vector<vector<FieldType>> ujd_transpose;
    transpose(ujd, ujd_transpose);

    for(auto it = ujd_transpose.begin(); it != ujd_transpose.end(); it++) {
        if(!checkConsistency(*it, degree)) {
            happiness = false;
            break;
        }
    }

    if(!isHappy(happiness)) { return false; }

    //recovered ui;
    vector<FieldType> ui;
    for(auto it = ujd_transpose.begin(); it != ujd_transpose.end(); it++) {
        ui.push_back(interpolate(*it, degree));
    }
    
    //broadcast computed ui's
    vector<byte> uiSend;
    vector<vector<byte>> uiRec(N1);
    encodeFieldElts(ui, uiSend);
    for(auto it = uiRec.begin(); it != uiRec.end(); it++) {
        it->resize(uiSend.size());
    }

    //uisend should be constant
    // = s1 + s2b + s3b^2 + ... where b = group elt and s1 is true secret. 
    //roundFunctionSyncBroadcast(uiSend, uiRec);
    for(int i = 0; i < inputs[1][partyID].size(); i++) {
        if(inputs[1][partyID][i] != uiSend) {
            return false;
        }
    }

    uiRec = outputs[1][partyID];


    /*
    vector<vector<FieldType>> uj(N1);

    for(int i = 0; i < activePartyIDs.size(); it++) {
        decodeFieldElts(uiRec[i], uj[i]);
    }

    vector<vector<FieldType>> uj_transpose;
    transpose(uj, uj_tranpose);
    
    //uiSend should be repeated field elt
    //uiRec should be 
    for(auto it = uj_transpose.begin(); it != uj_transpose.end(); it++) {
        if(!checkConsistency(*it, bigT - 1)) {
            happiness = false;
            break;
        }
    }

    if(!isHappy(happiness)) { return false; }

    //ui -> si values to coefficients
    for(auto it = uj_transpose.begin(); it != uj_transpose.end(); it++) {
        it->resize(bigT);
        vector<FieldType> alpha = getAlpha();
        vector<FieldType> si;
        lagrange(alpha, *it, si);
        for(auto it2 = si.begin(); it2 != si.end(); it2++) { reconstructedSecrets.push_back(*it2); }
    } */
    return true;
}

template <class FieldType>
pair<int,int> ProtocolParty<FieldType>::reconstructPublicVerifyAll(
    vector<FieldType>& secretShares, int degree, 
    vector<vector<vector<byte>>>& inputs, vector<vector<vector<byte>>>& outputs) {


    pair<int,int> result;

    vector<vector<vector<vector<byte>>>> inputsA(2);
    vector<vector<vector<vector<byte>>>> outputsA(2);


    vector<byte> secretShareBytes;

    encodeFieldElts(secretShares, secretShareBytes);

    vector<vector<byte>> secretShareBytesA;
    //is judge
    if(m_partyId == activePartyIDs[0]) {

        //exchange coefficients
        roundFunctionSyncForP1(secretShareBytes, secretShareBytesA);

        int wronground;

        for(wronground = 0; wronground < 2; wronground++) {
            result = JudgeIOConsistency(inputs[wronground], outputs[wronground], inputsA[wronground], outputsA[wronground]);
            if(result.first >= 0) { break; }
        }

        vector<vector<FieldType>> secretSharesA(N1);
        for(int i = 0; i < N1; i++) {
            decodeFieldElts(secretShareBytesA[i], secretSharesA[i]);
        }

        //TODO: check secretShares are consistent, if not, return self, inconsistent person. 
        for(int i = 0; i < secretShares.size(); i++) {
            vector<FieldType> secretI(N1);
            for(int j = 0; j < N1; j++) {
                secretI.push_back(secretSharesA[j][i]);
            }
            //check consistency of j; if non-consistent, return self, nonconsistent element
        }

        for(int i = 1; i < N1; i++) {
            if(result.first >= 0) { break; }

            if(!reconstructPublicVerifyOne(secretSharesA[i], degree, i, inputsA, outputsA)) 
            { //info provided by agent i conflicts with protocol
                result = make_pair(0, i);
            }
        }

        vector<byte> resPair = encodePair1(result);

        sendFromP1(resPair);

        //make some bounds checks?
        if(result.first > 0 && result.second > 0) { //pair does not include judge, others need to provide agree / disagree
            //broadcast all ij != ji messages
            vector<byte> comnum;
            comnum.push_back((byte)wronground);
            sendFromP1(comnum);
            sendFromP1(inputsA[wronground][result.first][result.second]);
            sendFromP1(outputsA[wronground][result.second][result.first]);

            vector<byte> disagree1(1);
            vector<byte> disagree2(1);

            //thie needs to be byzantine too.
            activeParties[result.first]->getChannel()->read(disagree1.data(), disagree1.size());
            activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());
            if(disagree1[0]) {
                removeParties(0, result.first);
            } else if(disagree2[0]){
                removeParties(0, result.second);
            } else {
                removeParties(result.first, result.second);
            }
        }

    } else { //conditions -1 = no inconsistency, if judge is inconsistent, then 
        activeParties[0]->getChannel()->write(secretShareBytes.data(), secretShareBytes.size());

        int wronground;

        for(wronground = 0; wronground < 2; wronground++) {
            result = NonJudgeIOConsistency(inputs[wronground], outputs[wronground]);
            if(result.first >= 0) { break; }
        }

        vector<byte> resultBytes(8);
        activeParties[0]->getChannel()->read(resultBytes.data(), resultBytes.size());
        result = decodePair1(resultBytes);
        
        if(result.first > 0 && result.second > 0) { 
            vector<byte> comnum(1);
            vector<byte> disagree1(1);
            vector<byte> disagree2(1);
            activeParties[0]->getChannel()->read(comnum.data(), comnum.size());
            vector<byte> msg(inputs[comnum[0]][result.first].size());
            
            if(activePartyIDs[result.first] == m_partyId) {
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                disagree1[0] = (msg != inputs[comnum[0]][result.second]);
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                sendFromP1(disagree1);
                activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());

            } else if(activePartyIDs[result.second] == m_partyId) {
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                disagree2[0] = (msg != outputs[comnum[0]][result.first]);
                activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());
                sendFromP1(disagree2);

            } else {
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                activeParties[result.first]->getChannel()->read(disagree1.data(), disagree1.size());
                activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());
            }

            if(disagree1[0]) {
                removeParties(0, result.first);
            } else if(disagree2[0]){
                removeParties(0, result.second);
            } else {
                removeParties(result.first, result.second);
            }
        }
    }

    return result;

}

template <class FieldType>
bool ProtocolParty<FieldType>::doubleShareRandom(int degree1, int degree2, vector<tuple<FieldType, FieldType>>& randomDoubles) {
    bool happiness = true;
    int eltSize = field->getElementSizeInBytes();
    vector<FieldType> d1share;
    vector<FieldType> d2share;
    vector<FieldType> d1coefficients;
    vector<FieldType> d2coefficients;
    FieldType secret = field->Random();
    d1coefficients.push_back(secret);
    d2coefficients.push_back(secret);
    for(int i = 0; i < degree1; i++) { d1coefficients.push_back(field->Random()); }
    for(int i = 0; i < degree2; i++) { d2coefficients.push_back(field->Random()); }

    vector<FieldType> alpha = getAlpha();
    
    for(int i = 0; i < alpha.size(); i++) {
        d1share.push_back(evaluatePolynomial(alpha[i], d1coefficients));
        d2share.push_back(evaluatePolynomial(alpha[i], d2coefficients));
    }

    vector<vector<byte>> d1sharebytes(N1);
    vector<vector<byte>> d2sharebytes(N1);
    vector<vector<byte>> d1shareallbytes(N1);
    vector<vector<byte>> d2shareallbytes(N1);
    
    for(int i = 0; i < N1; i++) {
        d1sharebytes[i].resize(eltSize);
        d2sharebytes[i].resize(eltSize);
        field->elementToBytes(d1sharebytes[i].data(),d1share[i]);
        field->elementToBytes(d2sharebytes[i].data(),d2share[i]);
        d1shareallbytes[i].resize(eltSize);
        d2shareallbytes[i].resize(eltSize);
    }

    roundFunctionSync(d1sharebytes, d1shareallbytes, 14);
    roundFunctionSync(d2sharebytes, d2shareallbytes, 15);


    vector<FieldType> d1shareall;
    vector<FieldType> d2shareall;

    for(int i = 0; i < N1; i++) {
        d1shareall.push_back(field->bytesToElement(d1shareallbytes[i].data()));
        d2shareall.push_back(field->bytesToElement(d2shareallbytes[i].data()));
    }

    vector<FieldType> r1shares(N1);
    vector<FieldType> r2shares(N1);

    HIM<FieldType> M;
    M.allocate(N1,N1,field);
    M.InitHIM();

    M.MatrixMult(d1shareall, r1shares);
    M.MatrixMult(d2shareall, r2shares);

    vector<vector<byte>> partialr1sharebytes(N1);
    vector<vector<byte>> partialr2sharebytes(N1);

    for(int i = 0; i < N1; i++) {
        partialr1sharebytes[i].resize(eltSize);
        partialr2sharebytes[i].resize(eltSize);
        if(i >= bigT) {
            field->elementToBytes(partialr1sharebytes[i].data(), r1shares[i]);
            field->elementToBytes(partialr2sharebytes[i].data(), r2shares[i]);
        } else {
            field->elementToBytes(partialr1sharebytes[i].data(), *(field->GetZero()));
            field->elementToBytes(partialr2sharebytes[i].data(), *(field->GetZero()));
        }
    }

    vector<vector<byte>> v1sharebytes(N1);
    vector<vector<byte>> v2sharebytes(N1);    
    for(int j = 0; j < activeParties.size(); j++) {
        v1sharebytes[j].resize(eltSize);
        v2sharebytes[j].resize(eltSize);
    }

    roundFunctionSync(partialr1sharebytes, v1sharebytes, 16);
    roundFunctionSync(partialr2sharebytes, v2sharebytes, 17);

    vector<FieldType> v1share(N1);
    vector<FieldType> v2share(N1);
    for(int j = 0; j < N1; j++) {
        v1share[j] = field->bytesToElement(v1sharebytes[j].data());
        v2share[j] = field->bytesToElement(v2sharebytes[j].data());
    }

    if(!checkConsistency(v1share, degree1) || !checkConsistency(v2share, degree2) || 
        interpolate(v1share, degree1) != interpolate(v2share, degree2)) {
        happiness = false;
    }

    //fault localization
    if(!isHappy(happiness)) { 

        vector<vector<vector<byte>>> inputs, outputs;
        inputs.push_back(d1sharebytes);
        inputs.push_back(d2sharebytes);
        inputs.push_back(partialr1sharebytes);
        inputs.push_back(partialr2sharebytes);
        outputs.push_back(d1shareallbytes);
        outputs.push_back(d2shareallbytes);
        outputs.push_back(v1sharebytes);
        outputs.push_back(v2sharebytes);

        doubleShareRandomVerifyAll(degree1, degree2, 
            d1coefficients, d2coefficients, inputs, outputs);

        return false; 

    }

    //if everyone is happy, output all remaining pairings.
    for(int i = 0; i < bigT; i++) {
        randomDoubles.push_back(make_tuple(r1shares[i], r2shares[i]));
    }
    return true;

}


template <class FieldType>
bool ProtocolParty<FieldType>::doubleShareRandomVerifyOne(int degree1, int degree2, 
    vector<FieldType>& d1coefficients, vector<FieldType>& d2coefficients, int partyID, 
    vector<vector<vector<vector<byte>>>>& inputs, vector<vector<vector<vector<byte>>>>& outputs) {

/*
    vector<vector<byte>>& d1sharebytes1, vector<vector<byte>>& d1sharebytesall1,
    vector<vector<byte>>& d2sharebytes1, vector<vector<byte>>& d2sharebytesall1,
    vector<vector<byte>>& partialr1sharebytes1, vector<vector<byte>>& v1sharebytes1,
    vector<vector<byte>>& partialr2sharebytes1, vector<vector<byte>>& v2sharebytes1) {
*/

    if(d1coefficients.size() != degree1 + 1 || d2coefficients.size() != degree2 + 1 
        || d1coefficients[0] != d2coefficients[0]) {
        return false;
    }

    vector<vector<byte>>& d1sharebytes1 = inputs[0][partyID];
    vector<vector<byte>>& d2sharebytes1 = inputs[1][partyID];
    vector<vector<byte>>& partialr1sharebytes1 = inputs[2][partyID];
    vector<vector<byte>>& partialr2sharebytes1 = inputs[3][partyID];

    vector<vector<byte>>& d1sharebytesall1 = outputs[0][partyID];
    vector<vector<byte>>& d2sharebytesall1 = outputs[1][partyID];
    vector<vector<byte>>& v1sharebytes1 = outputs[2][partyID];
    vector<vector<byte>>& v2sharebytes1 = outputs[3][partyID];

    int eltSize = field->getElementSizeInBytes();
    vector<FieldType> d1share;
    vector<FieldType> d2share;
    FieldType secret = field->Random();

    vector<FieldType> alpha = getAlpha();
    
    for(int i = 0; i < alpha.size(); i++) {
        d1share.push_back(evaluatePolynomial(alpha[i], d1coefficients));
        d2share.push_back(evaluatePolynomial(alpha[i], d2coefficients));
    }

    
    vector<vector<byte>> d1sharebytes(N1);
    vector<vector<byte>> d2sharebytes(N1);


    for(int i = 0; i < N1; i++) {
        d1sharebytes[i].resize(eltSize);
        d2sharebytes[i].resize(eltSize);
        field->elementToBytes(d1sharebytes[i].data(),d1share[i]);
        field->elementToBytes(d2sharebytes[i].data(),d2share[i]);
    }
    
    //generated wrong d1sharebytes
    if(d1sharebytes != d1sharebytes1 || d2sharebytes != d2sharebytes1) {
        return false;
    }

    vector<FieldType> d1shareall;
    vector<FieldType> d2shareall;

    for(int i = 0; i < N1; i++) {
        d1shareall[i] = field->bytesToElement(d1sharebytesall1[i].data());
        d2shareall[i] = field->bytesToElement(d2sharebytesall1[i].data());
    }

    vector<FieldType> r1shares(N1);
    vector<FieldType> r2shares(N1);

    HIM<FieldType> M;
    M.allocate(N1,N1,field);
    M.InitHIM();

    M.MatrixMult(d1shareall, r1shares);
    M.MatrixMult(d2shareall, r2shares);

    
    vector<vector<byte>> partialr1sharebytes(N1);
    vector<vector<byte>> partialr2sharebytes(N1);

    for(int i = 0; i < N1; i++) {
        partialr1sharebytes[i].resize(eltSize);
        partialr2sharebytes[i].resize(eltSize);
        if(i >= bigT) {
            field->elementToBytes(partialr1sharebytes[i].data(), r1shares[i]);
            field->elementToBytes(partialr2sharebytes[i].data(), r2shares[i]);
        } else {
            field->elementToBytes(partialr1sharebytes[i].data(), *(field->GetZero()));
            field->elementToBytes(partialr2sharebytes[i].data(), *(field->GetZero()));
        }
    }

    vector<vector<byte>> v1sharebytes(N1);
    vector<vector<byte>> v2sharebytes(N1);    
    for(int j = 0; j < N1; j++) {
        v1sharebytes[j].resize(eltSize);
        v2sharebytes[j].resize(eltSize);
    }

    //second set of communications wrong
    if(partialr1sharebytes != partialr1sharebytes1 || partialr2sharebytes != partialr2sharebytes1) {
        return false;
    }

    /*
    vector<FieldType> v1share(N1);
    vector<FieldType> v2share(N1);
    for(int j = 0; j < N1; j++) {
        v1share[j] = field->bytesToElement(v1sharebytes[j].data());
        v2share[j] = field->bytesToElement(v2sharebytes[j].data());
    }

    //we don't need to check all 0 for < bigT, since there are <= degree of them
    if(!checkConsistency(v1share, degree1) || !checkConsistency(v2share, degree2) || 
        interpolate(v1share, degree1) != interpolate(v2share, degree2)) {
        happiness = false;
    }

    if(!isHappy(happiness)) { return false; }
    */
    return true;
}

//called by everyone
template <class FieldType>
pair<int,int> ProtocolParty<FieldType>::doubleShareRandomVerifyAll(int degree1, int degree2, 
    vector<FieldType>& d1coefficients, vector<FieldType>& d2coefficients,
    vector<vector<vector<byte>>>& inputs, vector<vector<vector<byte>>>& outputs) {

    /*
    vector<vector<byte>>& d1sharebytes, vector<vector<byte>>& d1sharebytesall,
    vector<vector<byte>>& d2sharebytes, vector<vector<byte>>& d2sharebytesall,
    vector<vector<byte>>& partialr1sharebytes, vector<vector<byte>>& v1sharebytes,
    vector<vector<byte>>& partialr2sharebytes, vector<vector<byte>>& v2sharebytes) {
    */

    pair<int,int> result;

    vector<vector<vector<vector<byte>>>> inputsA(4);
    vector<vector<vector<vector<byte>>>> outputsA(4);


    /*
    //maybe make this into a single vector<vector<vector<vector<byte>>>>>?
    vector<vector<vector<byte>>> d1sharebytesA;
    vector<vector<vector<byte>>> d2sharebytesA;
    vector<vector<vector<byte>>> d1sharebytesallA;
    vector<vector<vector<byte>>> d2sharebytesallA;
    vector<vector<vector<byte>>> partialr1sharebytesA;
    vector<vector<vector<byte>>> partialr2sharebytesA;
    vector<vector<vector<byte>>> v1sharebytesA;
    vector<vector<vector<byte>>> v2sharebytesA;
    */

    vector<byte> d1coefficientBytes;
    vector<byte> d2coefficientBytes;

    encodeFieldElts(d1coefficients, d1coefficientBytes);
    encodeFieldElts(d2coefficients, d2coefficientBytes);

    vector<vector<byte>> d1coeffBytesA;
    vector<vector<byte>> d2coeffBytesA;

    //is judge
    if(m_partyId == activePartyIDs[0]) {

        //exchange coefficients
        roundFunctionSyncForP1(d1coefficientBytes, d1coeffBytesA);
        roundFunctionSyncForP1(d2coefficientBytes, d2coeffBytesA);

        int wronground;

        for(wronground = 0; wronground < 4; wronground++) {
            result = JudgeIOConsistency(inputs[wronground], outputs[wronground], inputsA[wronground], outputsA[wronground]);
            if(result.first >= 0) { break; }
        }

        for(int i = 1; i < N1; i++) {
            if(result.first >= 0) { break; }
            vector<FieldType> d1coefficientsI;
            vector<FieldType> d2coefficientsI;

            decodeFieldElts(d1coeffBytesA[i], d1coefficientsI);
            decodeFieldElts(d2coeffBytesA[i], d2coefficientsI); 

            if(!doubleShareRandomVerifyOne(degree1, degree2, d1coefficientsI, d2coefficientsI, i, inputsA, outputsA)) 
            { //info provided by agent i conflicts with protocol
                result = make_pair(0, i);
            }
        }

        vector<byte> resPair = encodePair1(result);

        sendFromP1(resPair);

        if(result.first > 0 && result.second > 0) { //pair does not include judge, others need to provide agree / disagree
            //broadcast all ij != ji messages

            vector<byte> comnum;
            comnum.push_back((byte)wronground);

            sendFromP1(comnum);
            sendFromP1(inputsA[wronground][result.first][result.second]);
            sendFromP1(outputsA[wronground][result.second][result.first]);

            vector<byte> disagree1(1);
            vector<byte> disagree2(1);

            //thie needs to be byzantine too.
            activeParties[result.first]->getChannel()->read(disagree1.data(), disagree1.size());
            activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());
            if(disagree1[0]) {
                removeParties(0, result.first);
            } else if(disagree2[0]){
                removeParties(0, result.second);
            } else {
                removeParties(result.first, result.second);
            }
        }

    } else { //conditions -1 = no inconsistency, if judge is inconsistent, then 
        activeParties[0]->getChannel()->write(d1coefficientBytes.data(), d1coefficientBytes.size());
        activeParties[0]->getChannel()->write(d2coefficientBytes.data(), d2coefficientBytes.size());

        int wronground;

        for(wronground = 0; wronground < 4; wronground++) {
            result = NonJudgeIOConsistency(inputs[wronground], outputs[wronground]);
            if(result.first >= 0) { break; }
        }

        vector<byte> resultBytes(8);
        activeParties[0]->getChannel()->read(resultBytes.data(), resultBytes.size());
        result = decodePair1(resultBytes);
        
        if(result.first > 0 && result.second > 0) { 
            vector<byte> comnum(1);
            vector<byte> disagree1(1);
            vector<byte> disagree2(1);
            activeParties[0]->getChannel()->read(comnum.data(), comnum.size());
            vector<byte> msg(inputs[comnum[0]][result.first].size());
            
            if(activePartyIDs[result.first] == m_partyId) {
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                disagree1[0] = (msg != inputs[comnum[0]][result.second]);
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                sendFromP1(disagree1);
                activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());

            } else if(activePartyIDs[result.second] == m_partyId) {
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                disagree2[0] = (msg != outputs[comnum[0]][result.first]);
                activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());
                sendFromP1(disagree2);

            } else {
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                activeParties[0]->getChannel()->read(msg.data(), msg.size());
                activeParties[result.first]->getChannel()->read(disagree1.data(), disagree1.size());
                activeParties[result.second]->getChannel()->read(disagree2.data(), disagree2.size());
            }

            if(disagree1[0]) {
                removeParties(0, result.first);
            } else if(disagree2[0]){
                removeParties(0, result.second);
            } else {
                removeParties(result.first, result.second);
            }
        }
    }

    return result;

}

//work refereeing into end of function they're in;
//will restart until succeed;
template <class FieldType>
void ProtocolParty<FieldType>::generateTriples(vector<tuple<FieldType, FieldType, FieldType>>& randomTriples) {
    while(selfActive) {
        vector<tuple<FieldType, FieldType>> a;
        vector<tuple<FieldType, FieldType>> b;
        vector<tuple<FieldType, FieldType>> r;
        vector<FieldType> d;
        vector<FieldType> dtrue;

        if( 
            !doubleShareRandom(T, T1, a) ||
            !doubleShareRandom(T, T1, b) ||
            !doubleShareRandom(T, 2*T1, r)
        ) {
            continue;
        }

        //cout << "Doubles Generated: " << T << " "<< T1 << endl;

        /*
        for(int i = 0; i < a.size(); i++) {
            cout << m_partyId << ": (" << get<0>(a[i]) << ", " << get<1>(a[i]) << ") (" << get<0>(b[i]) << ", " << get<1>(b[i]) << ") ("<< get<0>(r[i]) << ", " << get<1>(r[i]) << ")" << endl;
        }*/

        for(int i = 0; i < a.size(); i++) {
            d.push_back(get<1>(a[i]) * get<1>(b[i]) - get<1>(r[i]));
        }

        if(!reconstructPublic(d, 2 * T1, dtrue)) { continue; }

        for(int i = 0; i < a.size(); i++) {
            randomTriples.push_back(make_tuple(get<0>(a[i]), get<0>(b[i]), get<0>(r[i]) + dtrue[i]));
        }

        //tuple<FieldType, FieldType, FieldType> lastTriple = randomTriples[randomTriples.size() - 1];
        //cout << m_partyId << ": (" << get<0>(lastTriple) << ", " << get<1>(lastTriple) << ", "<< get<2>(lastTriple) <<")"<< endl;

        return;
    }

}

//!!! By Index, not id
template <class FieldType>
void ProtocolParty<FieldType>::removeParties(int party1, int party2) {
    N1 -= 2;
    T1 -= 1;
    if(party1 > party2) { swap(party1, party2); } //delete higher index first
    if(activePartyIDs[party1] == m_partyId || activePartyIDs[party2] == m_partyId) { selfActive = false; }
    activeParties.erase(activeParties.begin() + party2);
    activeParties.erase(activeParties.begin() + party1);
    activePartyIDs.erase(activePartyIDs.begin() + party2);
    activePartyIDs.erase(activePartyIDs.begin() + party1);

}

//checks if IO is consistent - if not, returns index pair of inconsistent communications. Otherwise, returns -1,-1
template <class FieldType>
pair<int,int> ProtocolParty<FieldType>::IOConsistency(vector<vector<vector<byte>>>& sentIO, vector<vector<vector<byte>>>& receivedIO) {
    //sentIO[i][j] is the vector<byte> sent by part i to party j
    //receiveIO[j][i] is the vector<byte> received by party j from party 
    for(int i = 0; i < N1; i++) {
        for(int j = 0; j < N1; j++) {
            if(i != j && activePartyIDs[i] != m_partyId && sentIO[i][j] != receivedIO[j][i]) {
                return make_pair(i,j);
            }
        }
    }
    return make_pair(-1,-1);
}

//if double, check degree is equal and and define d, d' shares;
//broadcast format - exchange 
//if degree2 < 0, then then no equality check needed.
//only call for judger, others call next fn
template <class FieldType>
pair<int,int> ProtocolParty<FieldType>::JudgeIOConsistency(vector<vector<byte>>& mySend, vector<vector<byte>>& myReceive, 
    vector<vector<vector<byte>>>& allSends, vector<vector<vector<byte>>>& allRecs) {

    if(activePartyIDs[0] != m_partyId) {
        cout << "Not Judge !!!" << endl;
        return make_pair(-1,-1);
    }

    int dim1 = mySend.size(); // == N1
    int dim2 = mySend[0].size(); // message size
    vector<byte> send1;
    vector<byte> rec1;

    vector<vector<byte>> allSends1;
    vector<vector<byte>> allRecs1;

    for(int i = 0; i < dim1; i++) {
        for(int j = 0; j < dim2; j++) {
            send1.push_back(mySend[i][j]);
            rec1.push_back(myReceive[i][j]);
        }
    }

    allSends.resize(N1);
    allRecs.resize(N1);
    for(int i = 0; i < N1; i++) {
        allSends1[i].resize(dim1 * dim2);
        allRecs1[i].resize(dim1 * dim2);
        allSends[i].resize(dim1);
        allRecs[i].resize(dim1);
    }

    roundFunctionSyncForP1(send1, allSends1);
    roundFunctionSyncForP1(rec1, allRecs1);

    //unflatten;
    for(int i = 0; i < N1; i++) {
        for(int j = 0; j < dim1; j++) {
            for(int k = 0; k < dim2; k++) {
                allSends[i][j].push_back(allSends1[i][dim2 * j + k]);
                allRecs[i][j].push_back(allRecs1[i][dim2 * j + k]);
            }
        }
    }

    pair<int, int> firstwrong = IOConsistency(allSends, allRecs);

    vector<byte> results = encodePair1(firstwrong);

    sendFromP1(results);

    return firstwrong;
}

//polynomial
template <class FieldType>
pair<int,int> ProtocolParty<FieldType>::NonJudgeIOConsistency(vector<vector<byte>>& mySend, vector<vector<byte>>& myReceive)  {

    if(activePartyIDs[0] == m_partyId) {
        cout << "Should be Judge !!!" << endl;
        return make_pair(-1,-1);
    }    

    int dim1 = mySend.size(); // == N1
    int dim2 = mySend[0].size(); // message size
    vector<byte> send1;
    vector<byte> rec1;

    for(int i = 0; i < dim1; i++) {
        for(int j = 0; j < dim2; j++) {
            send1.push_back(mySend[i][j]);
            rec1.push_back(myReceive[i][j]);
        }
    }

    int refId = activePartyIDs[0];
    activeParties[0]->getChannel()->write(send1.data(), send1.size());
    activeParties[0]->getChannel()->write(rec1.data(), rec1.size());

    vector<byte> results(8);

    //read from sent P1.
    activeParties[0]->getChannel()->read(results.data(), results.size());

    return decodePair1(results);
}

template <class FieldType>
int ProtocolParty<FieldType>::processMultiplications(int lastMultGate)
{

    int numProcessed = 0;
    vector<FieldType> a,b,c,d,e,dr,er;
    vector<int> outind;
    
    for (int k = circuit.getLayers()[currentCirciutLayer]; 
    k < circuit.getLayers()[currentCirciutLayer + 1]; k++)//go over only the logit gates
    {
        auto gate = circuit.getGates()[k];


        if (gate.gateType == MULT) {
            FieldType a1 = get<0>(randomTriplesArray[randomSharesOffset]);
            FieldType b1 = get<1>(randomTriplesArray[randomSharesOffset]);
            FieldType c1 = get<2>(randomTriplesArray[randomSharesOffset]);

            a.push_back(a1);
            b.push_back(b1);
            c.push_back(c1);
            d.push_back(gateShareArr[gate.input1] - a1);
            e.push_back(gateShareArr[gate.input2] - b1);
            outind.push_back(gate.output);

            numProcessed++;

        }
    }

    //TODO: this should be robust rather than restartable (or we can just generate more triples)
    reconstructPublic(d, T, dr);
    reconstructPublic(e, T, er);

    for(int i = 0; i < numProcessed; i++) {
        gateShareArr[outind[i]] = dr[i] * er[i] + dr[i] * b[i] + er[i] * a[i] + c[i];
        cout << "Party: " << m_partyId << " Gate: m" << outind[i] << " output: " << gateShareArr[outind[i]] << endl;
    }
    return numProcessed;
}


template <class FieldType>
int ProtocolParty<FieldType>::processNotMult(){
    int count=0;
    cout << m_partyId << "gateShareArr1";
    printVector(gateShareArr);
    for(int k=circuit.getLayers()[currentCirciutLayer]; k < circuit.getLayers()[currentCirciutLayer+1]; k++)
    {

        string gateType = "Unknown";
        // add gate
        if(circuit.getGates()[k].gateType == ADD)
        {
            gateType = "Add";
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] + gateShareArr[circuit.getGates()[k].input2];
            count++;
        }

        else if(circuit.getGates()[k].gateType == SUB)//sub gate
        {
            gateType = "Sub";
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] - gateShareArr[circuit.getGates()[k].input2];
            count++;
        }
        else if(circuit.getGates()[k].gateType == SCALAR)
        {
            gateType = "Scalar";
            long scalar(circuit.getGates()[k].input2);
            FieldType e = field->GetElement(scalar);
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] * e;
            count++;
        }
        else if(circuit.getGates()[k].gateType == SCALAR_ADD)
        {
            gateType = "ScalarAdd";
            long scalar(circuit.getGates()[k].input2);
            FieldType e = field->GetElement(scalar);
            gateShareArr[circuit.getGates()[k].output] = gateShareArr[circuit.getGates()[k].input1] + e;
            count++;
        }
        else if(circuit.getGates()[k].gateType == RANDOM) 
        {
            gateType = "Random";
            gateShareArr[circuit.getGates()[k].output] = get<0>(randomTriplesArray[randomSharesOffset]);
            randomSharesOffset++;
            count++;
        }

        cout << "Party: " << m_partyId << " Gate: " << circuit.getGates()[k].output << "index: " << k << " output: " << gateShareArr[circuit.getGates()[k].output] << gateType << endl;

    }

    cout << m_partyId << "gateShareArr2";
    printVector(gateShareArr);
    return count;

}

template <class FieldType>
void ProtocolParty<FieldType>::run() {

    for (iteration=0; iteration<times; iteration++){

        auto t1start = high_resolution_clock::now();
        timer->startSubTask("Offline", iteration);
        cout << "Offline" << endl;
        runOffline();
        timer->endSubTask("Offline", iteration);
        timer->startSubTask("Online", iteration);
        cout << "Online" << endl;
        runOnline();
        timer->endSubTask("Online", iteration);

        auto t2end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2end-t1start).count();
        protocolTimer->totalTimeArr[iteration] = duration;

        cout << "time in milliseconds for protocol: " << duration << endl;
    }


}

template <class FieldType>
void ProtocolParty<FieldType>::runOffline() {
    auto t1 = high_resolution_clock::now();
    timer->startSubTask("preparationPhase", iteration);
    if(preparationPhase() == false) {
        if(flag_print) {
            cout << "cheating!!!" << '\n';}
        return;
    }
    else {
        if(flag_print) {
            cout << "no cheating!!!" << '\n' << "finish Preparation Phase" << '\n';}
    }
    timer->endSubTask("preparationPhase", iteration);
    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds preparationPhase: " << duration << endl;
    }
    protocolTimer->preparationPhaseArr[iteration] =duration;
}

template <class FieldType>
void ProtocolParty<FieldType>::runOnline() {

    auto t1 = high_resolution_clock::now();
    timer->startSubTask("inputPhase", iteration);
    inputPhase();
    timer->endSubTask("inputPhase", iteration);
    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    protocolTimer->inputPreparationArr[iteration] = duration;
    if(flag_print_timings) {
        cout << "time in milliseconds inputPhase: " << duration << endl;
    }


    t1 = high_resolution_clock::now();
    timer->startSubTask("ComputePhase", iteration);
    computationPhase();
    timer->endSubTask("ComputePhase", iteration);
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    protocolTimer->computationPhaseArr[iteration] = duration;



    if(flag_print_timings) {
        cout << "time in milliseconds computationPhase: " << duration << endl;
    }

    t1 = high_resolution_clock::now();
    timer->startSubTask("outputPhase", iteration);
    outputPhase();
    timer->endSubTask("outputPhase", iteration);
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    protocolTimer->outputPhaseArr[iteration] = duration;

    if(flag_print_timings) {
        cout << "time in milliseconds outputPhase: " << duration << endl;
    }

}

template <class FieldType>
bool ProtocolParty<FieldType>::preparationPhase()
{

    randomSharesOffset = 0; //use this for triples instead
    int triplesNeeded = circuit.getNrOfInputGates() + 
    circuit.getNrOfMultiplicationGates() + circuit.getNrOfRandomGates();
    for(int gS = 0; gS < triplesNeeded; gS += bigT) {
        generateTriples(randomTriplesArray);
        
        if(flag_print_timings && m_partyId == activePartyIDs[0]) {
            if(((gS + bigT) % 20) < (gS % 20)) { cout << "Triples Generated: " << gS+bigT << "/" << triplesNeeded << endl; }
        }
    }

    return true;
}

template <class FieldType>
void ProtocolParty<FieldType>::initializationPhase()
{
    gateShareArr.resize((M - circuit.getNrOfOutputGates())); // my share of the gate (for all gates)
}

//input s, share r, reconstruct r towards user, then user broadcasts s-r, so s-r + r = s
template <class FieldType>
void ProtocolParty<FieldType>::inputPhase()
{
    int k = 0;
    int ind = 0;
    int index = 0;

    // prepare the shares for the input
    //TODO: batch this
    cout << "gateShareArr1";
    printVector(gateShareArr);

    while(k < numOfInputGates)
    {
        if(circuit.getGates()[ind].gateType == INPUT) {

            int inputParty = circuit.getGates()[ind].party;
            FieldType randomShare = get<0>(randomTriplesArray[randomSharesOffset]); //[r]
            randomSharesOffset++;
            FieldType output = reconstructPrivate(randomShare, T, inputParty); //r or zero if not m_partyId
            if (inputParty == m_partyId) {
                auto input = myInputs[index]; //maybe this gets messed up with ind or k? 
                index++;

                //need to make sure random offsets sync up between active / inactive parties
                FieldType s = field->GetElement(input); //s
                vector<byte> sminusrbytes(field->getElementSizeInBytes());
                FieldType smo = s - output;
                field->elementToBytes(sminusrbytes.data(), smo);

                sendFromP1(sminusrbytes);

                int gateshareindex = circuit.getGates()[ind].output;
                gateShareArr[gateshareindex] = randomShare + s - output;
                cout << "Party: " << m_partyId << " Gate: " << ind << " Index " << gateshareindex << " RandomInput: " << gateShareArr[circuit.getGates()[ind].output] << endl;
            } else if(selfActive) {

                vector<byte> sminusrbytes(field->getElementSizeInBytes());

                activeParties[inputParty]->getChannel()->read(sminusrbytes.data(), sminusrbytes.size());

                int gateshareindex = circuit.getGates()[ind].output;
                gateShareArr[gateshareindex] = randomShare + field->bytesToElement(sminusrbytes.data());
                cout << "Party: " << m_partyId << " Gate: " << ind << " Index " << gateshareindex << " RandomInput: " << gateShareArr[circuit.getGates()[ind].output] << endl;
            }
            k++;
            if(flag_print_timings) {
                { cout << "Inputs: " << k << "/" << numOfInputGates << endl; }
            }
        }
        
        ind++;
        if(ind >= M) {
            cout << "Ind too high! " << k << endl;
        }
    }

    cout << "gateShareArr";
    printVector(gateShareArr);

    size_t numgates = circuit.getGates().size();
    for(size_t i = 0; i < numgates; i++) {
        if(circuit.getGates()[i].gateType == INPUT) {
            cout << "Party: " << m_partyId << " Gate: i" << i << " output: " << gateShareArr[circuit.getGates()[i].output] << endl;
        }
    }
}

template <class FieldType>
void ProtocolParty<FieldType>::computationPhase() {
    int count = 0;
    int countNumMult = 0;
    int countNumMultForThisLayer = 0;

    int numOfLayers = circuit.getLayers().size();
    for(int i=0; i<numOfLayers-1;i++){
        if(flag_print) {
            cout << "Processing Layer: " << i << endl;
        }
        currentCirciutLayer = i;
        count = processNotMult();

        countNumMultForThisLayer = processMultiplications(countNumMult);//send the index of the current mult gate
        countNumMult += countNumMultForThisLayer;;
        count+=countNumMultForThisLayer;

    }
}

//TODO: make changes for activeParty system
template <class FieldType>
void ProtocolParty<FieldType>::outputPhase()
{
    int count=0;
    vector<FieldType> x1(N); // vector for the shares of my outputs
    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);
    vector<vector<byte>> recBufBytes(N);

    FieldType num;
    ofstream myfile;
    myfile.open(outputFile);

    for(int k=M-numOfOutputGates; k < M; k++)
    {
        if(circuit.getGates()[k].gateType == OUTPUT)
        {
            // send to party (which need this gate) your share for this gate
            sendBufsElements[circuit.getGates()[k].party].push_back(gateShareArr[circuit.getGates()[k].input1]);
        }
    }


    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);
        recBufBytes[i].resize(sendBufsElements[m_partyId].size()*fieldByteSize);
//        for(int j=0; j<sendBufsElements[i].size();j++) {
//            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
//        }


        encodeFieldElts(sendBufsElements[i], sendBufsBytes[i]);
    }



    //comm->roundfunctionI(sendBufsBytes, recBufBytes,7);
    roundFunctionSync(sendBufsBytes, recBufBytes,7);



    int counter = 0;
    if(flag_print) {
        cout << "endnend" << endl;}
    for(int k=M-numOfOutputGates ; k < M; k++) {
        if(circuit.getGates()[k].gateType == OUTPUT && circuit.getGates()[k].party == m_partyId)
        {
            for(int i=0; i < N; i++) {

                x1[i] = field->bytesToElement(recBufBytes[i].data() + (counter*fieldByteSize));
            }

            // my output: reconstruct received shares
            if (!checkConsistency(x1, T))
            {
                // someone cheated!
                //if(flag_print) {
                    cout << "cheating!!!" << '\n';//}
                return;
            }
            if(flag_print_output)
                cout << "the result for "<< circuit.getGates()[k].input1 << " is : " << field->elementToString(interpolate(x1, N1)) << '\n';


            counter++;
        }
    }

    // close output file
    myfile.close();
}

template <class FieldType>
void ProtocolParty<FieldType>::readMyInputs()
{

    //cout<<"inputs file" << inputsFile<<endl;
    ifstream myfile;
    long input;
    int i =0;
    myfile.open(inputsFile);
    do {
        myfile >> input;
        myInputs[i] = input;
        i++;
    } while(!(myfile.eof()));
    myfile.close();


}

template <class FieldType>
ProtocolParty<FieldType>::~ProtocolParty()
{
    protocolTimer->writeToFile();
    delete protocolTimer;
    delete field;
    delete timer;
    //delete comm;
}


#endif /* PROTOCOLPARTY_H_ */
