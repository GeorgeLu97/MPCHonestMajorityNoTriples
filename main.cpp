
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


int main(int argc, char* argv[])
{

    CmdParser parser;
    auto parameters = parser.parseArguments("", argc, argv);
    int times = stoi(parser.getValueByKey(parameters, "internalIterationsNumber"));


    string fieldType = parser.getValueByKey(parameters, "fieldType");
    cout<<"fieldType = "<<fieldType<<endl;

    int polynomial[] = {4, 1, 2, 3, 5}; // size = 5
    int alpha[] = {1, 4, 9, 20, 5, 8, 6, 3, 7, 2}; // size = 10
    
    vector<int> result;
    vector<int> poly(polynomial, polynomial+5);
    vector<int> alph(alpha, alpha+10);

    if(fieldType.compare("ZpMersenne31") == 0)
      {
        // testReconstruct<ZpMersenneIntElement>(poly, alph, result);
        ProtocolParty<ZpMersenneIntElement> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
            protocol.run();

        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

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

    return 0;
}
