Prime Numbers Factorization in JavaScript
=========================================

Overview
This project provides an efficient algorithm for prime factorization using trial division, Miller-Rabin primality testing, and Pollard's Rho algorithm. It can handle both small and large numbers effectively.

Features:
- Supports BigInt for large numbers.
- Uses Miller-Rabin for probabilistic primality testing.
- Implements Pollard’s Rho for efficient factorization.
- Optimized trial division for small factors.

Installation:
1. Clone the repository:
   git clone https://github.com/dixy52-beep/Prime-Numbers-Factorization-JavaScript.git
2. Navigate to the project folder:
   cd Prime-Numbers-Factorization-JavaScript

Usage:
- Run the script using Node.js:
  node script.js
- Modify `main()` to test different numbers.

Example:
console.log(primeFactorization(13195n));
// Output: [5n, 7n, 13n, 29n]

Contributions:
You are welcome to fork the repository and submit pull requests for enhancements or bug fixes.

License:
This project is licensed under the MIT License. See the LICENSE file for details.
