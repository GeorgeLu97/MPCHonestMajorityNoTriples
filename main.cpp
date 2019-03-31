#include <stdlib.h>
#include "ProtocolParty.h"
#include "BAParty.h"
#include "ECC.h"
#include "ZpKaratsubaElement.h"
#include <smmintrin.h>
#include <inttypes.h>
#include <stdio.h>
#include <x86intrin.h>


template <class FieldType>
void testReconstruct(const vector<int>& poly,
                     const vector<int>& alpha,
                     vector<int>& re_poly){

  ECC<FieldType> ecc;
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
    fieldY[i] = ecc.evalPolynomial( fieldX[i], fieldPoly );
    
  }
  
  ecc.setAlpha(fieldX);
  
  // output to check
  cout << "evaluation result:" << endl;
  for(int i=0; i<nPoints; i++){
    cout << "(" << fieldX[i] << " " << fieldY[i] << ") ";
  }
  cout << endl;

  // corrupt (11 - 5) / 2 = 3 points
  fieldY[1] = FieldType(0) - fieldY[2];
  fieldY[10] = FieldType(1) - fieldY[3];
  fieldY[5] = fieldY[0] - fieldY[4];

  // try reconstruction
  ecc.reconstruct(fieldY, nCoeff-1, result);
  
  cout << "evaluation after reconstruction:" << endl;
  for(int i=0; i<nPoints; i++){
    cout << "(" << fieldX[i] << " "
         << ecc.evalPolynomial( fieldX[i], result ) << ") ";
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
    /*
        // testReconstruct<ZpMersenneIntElement>(poly, alph, result);
      }
    else if(fieldType.compare("ZpMersenne61") == 0)
      {
        // testReconstruct<ZpMersenneLongElement>(poly, alph, result);
        testReconstruct<ZpMersenneLongElement>(poly, alph, result);
      }
    else if(fieldType.compare("ZpKaratsuba") == 0)
      {
        // testReconstruct<ZpKaratsubaElement>(poly, alph, result);
        testReconstruct<ZpKaratsubaElement>(poly, alph, result);
      }
    else if(fieldType.compare("GF2m") == 0)
      {
        // testReconstruct<GF2E>(poly, alph, result);
        testReconstruct<GF2E>(poly, alph, result);
      }
    else if(fieldType.compare("Zp") == 0)
      {
        // testReconstruct<ZZ_p>(poly, alph, result);
        testReconstruct<ZZ_p>(poly, alph, result);
      }
	*/

    return 0;
}
