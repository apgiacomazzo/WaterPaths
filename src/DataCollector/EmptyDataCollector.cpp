//
// Created by bernardoct on 8/26/17.
//

#include "EmptyDataCollector.h"
#include "../Utils/Constants.h"

EmptyDataCollector::EmptyDataCollector() : DataCollector
                                                   (Constants::NON_INITIALIZED,
                                                    nullptr,
                                                    Constants::NON_INITIALIZED,
                                                    0) {}

string EmptyDataCollector::printTabularString(int week) {
    string output = "";
    return output;
}

string EmptyDataCollector::printCompactString(int week) {
    string output = "";
    return output;
}

string EmptyDataCollector::printTabularStringHeaderLine1() {
    string output = "";
    return output;
}

string EmptyDataCollector::printTabularStringHeaderLine2() {
    string output = "";
    return output;
}

string EmptyDataCollector::printCompactStringHeader() {
    return nullptr;
}

void EmptyDataCollector::collect_data() {

}
