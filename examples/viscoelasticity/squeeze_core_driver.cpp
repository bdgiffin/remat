#include "arithmetic.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: squeeze_core_driver p q\n";
    return EXIT_FAILURE;
  }

  LongInteger p = 0;
  LongInteger q = 0;
  try {
    p = std::stoll(argv[1]);
    q = std::stoll(argv[2]);
  } catch (const std::exception &error) {
    std::cerr << "invalid p or q: " << error.what() << "\n";
    return EXIT_FAILURE;
  }
  if (p <= 0 || q <= 0) {
    std::cerr << "p and q must be positive\n";
    return EXIT_FAILURE;
  }

  LongInteger first = 0;
  LongInteger second = 0;
  while (std::cin >> first >> second) {
    auto [mapped_first, mapped_second] = squeezeE(first, second, p, q);
    std::cout << mapped_first << ' ' << mapped_second << '\n';
  }

  if (!std::cin.eof()) {
    std::cerr << "failed to parse input point stream\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
