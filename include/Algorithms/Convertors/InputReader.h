#ifndef INPUTREADER_H
#define INPUTREADER_H

#include <string>
#include <mlpack/core.hpp>
class InputReader
{
    public:
        InputReader();
        static void XYZtoMAT(std::string fl, arma::mat & data);

        //virtual ~InputReader();

    protected:

    private:
};

#endif // INPUTREADER_H
