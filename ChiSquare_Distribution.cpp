/*
  File:         ChiSquare_Distribution.cpp
  Version:      0.0.2
  Date:         27-Jan-2019
  Revision:     31-Jan-2019
  Author:       Jerome Drouin (jerome.p.drouin@gmail.com)

  Editions:	Please go to ChiSquare_Distribution.h for Edition Notes.

  ChiSquare_Distribution.cpp - Library for 'duino
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

*/


#include "Arduino.h"
#include "ChiSquare_Distribution.h"
#include <Gamma_Function.h>


// Constructor /////////////////////////////////////////////////////////////////
// Function that handles the creation and setup of instances
ChiSquare_Distribution::ChiSquare_Distribution(int _Nu)
{

	// Object parameter's error handling
	error = 1;
	
	if (_Nu<=0) 	error =	-2;		// Nu must be Positive
	if (_Nu%1 !=0) 	error = -3;		// Nu must be Integer

	//Set initial values	
	Nu			= _Nu;		// 

}


// Public Methods //////////////////////////////////////////////////////////////
//Probability Density Function
double ChiSquare_Distribution::GetPDF(double x)
{
	//if error, return error
	if (error<0) return error;

	if (x > 0) {
		return 1.0/(pow(2.0,0.5*Nu) * Gamma_Function((double)(Nu*.5))) * pow(x,0.5*Nu-1) * exp(-0.5*x);
	} else {
		return 0;
	}

}


//Cumulative Distribution Function
double ChiSquare_Distribution::GetCDF(double x)
{
	//if error, return error
	if (error<0) return error;

	if (x > 0) {
		return 1.0/Gamma_Function((double)(Nu*.5)) * Lower_Incomplete_Gamma_Function(x*0.5, (double)(Nu*.5));
	} else {
		return 0;
	}
}


//Mean
double ChiSquare_Distribution::GetMean(void)
{
	return Nu;
}


//Variance
double ChiSquare_Distribution::GetVariance(void)
{
	return 2.0*Nu;
}


//Std Deviation
double ChiSquare_Distribution::GetStdDeviation(void)
{
	return sqrt(2.0*Nu);
}


//Skewness
double ChiSquare_Distribution::GetSkewness(void)
{
	return sqrt(8.0/Nu);
}


//Kurtosis
double ChiSquare_Distribution::GetKurtosis(void)
{
	return 12.0/Nu;
}


//Return Quantile z(P) from probability P
double ChiSquare_Distribution::GetQuantile(double p)
{
double Vm;
double Vh = 100;
double Vl = 0;
double Pr;
int i = 0;
double Eps;

	//if error, return error
	if (error<0) return error;

	if (p <= 0.0) {
		return Vl;
	} else if (p >= 1.0) {
        	return Vh;
	} else {        
        	do 
		{
          		i++;
          		Vm = (Vh + Vl) / 2;
            
			Pr = GetCDF(Vm);
          		Eps = abs(Pr - p);
			        
          		//New Boundary selection
          		if (Pr > p) {
				Vh = Vm;
          		} else {
				Vl = Vm;
			}
            
        	} 
		while ((Eps > CONSTANT_EpsStop) && (i < 70));
	}
            
        if (i >= 70) {
            return -9999;
        } else {
            return Vm;
    	}

}


double ChiSquare_Distribution::GetNu(void)
{
	return Nu;
}


// Private Methods /////////////////////////////////////////////////////////////
// Functions only available to other functions in this library


// /////////////////////////////////////////////////////////////////////////////

