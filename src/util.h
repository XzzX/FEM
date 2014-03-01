#ifndef UTIL_H
#define UTIL_H

#include	<sstream>

/**
Converts a string to an arbitrary type. >> operator must be defined for the target type.
@param string string which should be converted
@return converted string
**/
template<typename T>
T StringTo(const std::string& string){
    T valor;

    std::stringstream stream(string);
    stream >> valor;
    if (stream.fail()) {
        std::runtime_error e(string);
        throw e;
    }
    return valor;
}

#endif