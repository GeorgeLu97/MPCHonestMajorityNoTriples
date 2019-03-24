
#include <stdlib.h>
#include "ProtocolParty.h"
#include "ZpKaratsubaElement.h"
#include <smmintrin.h>
#include <inttypes.h>
#include <stdio.h>
#include <x86intrin.h>
#include "BAParty.h"


template <class FieldType>
void testReconstruct(const vector<int>& poly,
                     const vector<int>& alpha,
                     vector<int>& re_poly){
  int nCoeff = poly.size();
  int nPoints = alpha.size();
  vector<FieldType> fieldPoly(nCoeff, FieldType());
  vector<FieldType> fieldX(nPoints, FieldType());
  vector<FieldType> fieldY(nPoints, FieldType());
  vector<FieldType> result;

  // build polynomial
  for(int i=0; i<nCoeff; i++){
    fieldPoly[i]= FieldType( poly[i] );
  }

  // build sample points
  for(int i=0; i<nPoints; i++){
    fieldX[i] = FieldType( alpha[i] );
    fieldY[i] = evalPolynomial<FieldType>( fieldX[i], fieldPoly );
    
  }

  // output to check
  cout << "evaluation result:" << endl;
  for(int i=0; i<nPoints; i++){
    cout << "(" << fieldX[i] << " " << fieldY[i] << ") ";
  }
  cout << endl;

  // try reconstruction
  reconstruct <FieldType> (fieldX, fieldY, nCoeff-1, result);
  
  cout << "evaluation after reconstruction:" << endl;
  for(int i=0; i<nPoints; i++){
    cout << "(" << fieldX[i] << " "
         << evalPolynomial<FieldType>( fieldX[i], result ) << ") ";
  }
  cout << endl;
  
  return;
}

// a BA protocol as a subprotocol

int main(int argc, char* argv[])
{

    CmdParser parser;
    auto parameters = parser.parseArguments("", argc, argv);
    int times = stoi(parser.getValueByKey(parameters, "internalIterationsNumber"));


    string fieldType = parser.getValueByKey(parameters, "fieldType");
    cout<<"fieldType = "<<fieldType<<endl;

    if(fieldType.compare("ZpMersenne31") == 0)
    {
        ProtocolParty<ZpMersenneIntElement> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
            protocol.run();

        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

    }
    else if(fieldType.compare("ZpMersenne61") == 0)
    {

        ProtocolParty<ZpMersenneLongElement> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
        protocol.run();
        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

    }

    else if(fieldType.compare("ZpKaratsuba") == 0) {
        ProtocolParty<ZpKaratsubaElement> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
        protocol.run();

        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2 - t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';
    }



    else if(fieldType.compare("GF2m") == 0)
    {
        ProtocolParty<GF2E> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
            protocol.run();

        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';
    }

    else if(fieldType.compare("Zp") == 0)
    {
    	
    	cout << "FieldType Zp" << endl;
        ProtocolParty<ZZ_p> protocol(argc, argv);

        cout << "ProtocolParty Init" << endl;

        auto t1 = high_resolution_clock::now();

        protocol.run();

        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

    }

    return 0;
}
