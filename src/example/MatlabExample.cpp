#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include <iostream>
void callSQRT() {

    using namespace matlab::engine;

    // connect to MATLAB engine synchronously
    std::unique_ptr<MATLABEngine> matlabPtr = connectMATLAB();   

    //Create MATLAB data array factory
    matlab::data::ArrayFactory factory;

    // Define a four-element typed array
    matlab::data::TypedArray<double> const argArray = 
        factory.createArray({ 1,4 }, { -2.0, 2.0, 6.0, 8.0 });

    // Call MATLAB sqrt function on the data array
    matlab::data::Array const results = matlabPtr->feval(u"sqrt", argArray);

    // Display results
    for (int i = 0; i < results.getNumberOfElements(); i++) {
        double a = argArray[i];
        std::complex<double> v = results[i];
        double realPart = v.real();
        double imgPart = v.imag();
        std::cout << "Square root of " << a << " is " << 
            realPart << " + " << imgPart << "i" << std::endl;
    }

    std::cout << "Test running script" << std::endl;
    matlabPtr->eval(u"run(\'/home/cm58349/Documents/OrbitalDetermination/src/example/MatRunTest.m\')");
}

int main() {
    callSQRT();
    return 0;
}