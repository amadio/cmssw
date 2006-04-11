//   COCOA class implementation file
//Id:  ALIRmDataFromFile.cc
//CAT: Model
//
//   History: v1.0 
//   Pedro Arce

#include "Alignment/CocoaModel/interface/ALIRmDataFromFile.h"
#include "Alignment/CocoaUtilities/interface/ALIUtils.h"


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ALIRmDataFromFile::ALIRmDataFromFile()
{
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ALIbool ALIRmDataFromFile::setAngle( const ALIstring& coord, const ALIdouble val ){

  if( coord == "X" ) {
    return setAngleX( val );
  }else if( coord == "Y" ) {
    return setAngleY( val );
  }else if( coord == "Z" ) {
    return setAngleZ( val );
  }else {   
    std::cerr << "!!! FATAL ERROR ALIRmDataFromFile::setAngle. Coordinate must be X, Y or Z, it ii " << coord << std::endl;
    std::exception();
  }
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ALIbool ALIRmDataFromFile::setAngleX( const ALIdouble val )
{
  theAngleX = val;
  theDataFilled += "X";
  return 1;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ALIbool ALIRmDataFromFile::setAngleY( const ALIdouble val )
{
  theAngleY = val;
  theDataFilled += "Y";
  return 1;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ALIbool ALIRmDataFromFile::setAngleZ( const ALIdouble val )
{
  theAngleZ = val;
  theDataFilled += "Z";
  return 1;
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void ALIRmDataFromFile::constructRm()
{
  if( theDataFilled.find("X") == -1 ||  theDataFilled.find("Y") == -1 ||  theDataFilled.find("Z") == -1 ){
    std::cerr << "!!!  ALIRmDataFromFile::constructRm. FATAL ERROR: building rm while one angle is missing: " << theDataFilled << std::endl;
  } else {
    theRm = CLHEP::HepRotation();
    theRm.rotateX( theAngleX );
    theRm.rotateY( theAngleY );
    theRm.rotateZ( theAngleZ );
  }

}
