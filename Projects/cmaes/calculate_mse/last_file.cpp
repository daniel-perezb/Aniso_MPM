#include "functions.h"
#include <fstream>

// Function to find the last file number instead of hyperparameter
int findLastFileNumber() {
  int lastFileNumber = 0;
  for (int i = 0; i < 1000; i++) {
    std::string filename =
        "output/RIG_fleece_fibres/data_" + std::to_string(i) + ".dat";
    std::ifstream file(filename);
    if (!file) {
      lastFileNumber = i - 1;
      break;
    }
  }
  return lastFileNumber;
}
