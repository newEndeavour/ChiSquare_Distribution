/*
  File:         ChiSquare_Distribution.h
  Version:      0.0.2
  Date:         27-Jan-2019
  Revision:     31-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  ChiSquare_Distribution.h - Library for 'duino
  https://github.com/newEndeavour/ChiSquare_Distribution

  ChiSquare_Distribution implements the ChiSquare distribution 
  with v (Nu) degrees of Freedom.

  Copyright (c) 2018-2019 Jerome Drouin  All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Editions:
  - 0.0.1	: First version
  - 0.0.2	: Bug corrections. Replacement of double variables for double

*/


// ensure this library description is only included once
#ifndef ChiSquare_Distribution_h
#define ChiSquare_Distribution_h


#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif


//ChiSquare Constant Parameters
#define CONSTANT_aa 0.2316419
#define CONSTANT_b1 1.330274429#define CONSTANT_b2 1.821255978
#define CONSTANT_b3 1.781477937
#define CONSTANT_b4 0.356563782
#define CONSTANT_b5 0.31938153
#define CONSTANT_EpsStop 0.0000001
    
//Number Pi
#define CONSTANT_Pi 3.14159265358979


// library interface description
class ChiSquare_Distribution
{
  // user-accessible "public" interface
  public:
  // methods
	ChiSquare_Distribution(int _Nu);
	
	double 	GetPDF(double x);
	double 	GetCDF(double x);
	double 	GetQuantile(double p);

	double 	GetNu(void);

  // library-accessible "private" interface
  private:
  // variables
	int 	error;

	int	Nu;			// 
	
  // methods
};

#endif
