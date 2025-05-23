// Helper function for BigInt modular exponentiation (a^b % m)
function bigIntPower(base, exp, mod) {
    let res = 1n;
    base %= mod;
    while (exp > 0n) {
        if (exp % 2n === 1n) res = (res * base) % mod;
        base = (base * base) % mod;
        exp /= 2n;
    }
    return res;
}

// Helper function for BigInt Greatest Common Divisor (GCD) using Euclidean algorithm
function bigIntGCD(a, b) {
    a = a < 0n ? -a : a;
    b = b < 0n ? -b : b;
    while (b) {
        [a, b] = [b, a % b];
    }
    return a;
}

// Helper function to generate a random BigInt in the range [min, max]
// This is a simple implementation using Math.random() and might have some bias.
// Sufficient for generating random bases in Miller-Rabin for probabilistic testing.
function randomBigInt(min, max) {
    if (min > max) {
        throw new Error("Invalid range for randomBigInt: min > max");
    }
    const range = max - min + 1n;
    if (range <= 0n) {
         return min;
    }

    // Estimate the number of bits needed for the range
    const rangeBits = range.toString(2).length;
    let randomVal = 0n;
    // Use a chunk size less than 53 to stay within Math.random's integer precision
    const chunkSize = 48;
    const numChunks = Math.ceil(rangeBits / chunkSize);

    // Generate random bits chunk by chunk
    for (let i = 0; i < numChunks; i++) {
        const chunk = BigInt(Math.floor(Math.random() * (2 ** chunkSize)));
        randomVal = (randomVal << BigInt(chunkSize)) | chunk;
    }

    // Reduce the random number to the desired range [0, range-1] and add min
    // This modulo approach can introduce bias, but is acceptable for MR bases.
    return (randomVal % range) + min;
}

// Miller-Rabin primality test for a single base 'a'.
// Returns false if n is composite, true if n is a strong pseudoprime to base 'a'.
// Assumes n is odd and n > 3.
function millerRabinTest(n, a) {
    let d = n - 1n;
    let s = 0n;
    while (d % 2n === 0n) {
        d /= 2n;
        s++;
    }

    let x = bigIntPower(a, d, n);

    if (x === 1n || x === n - 1n) return true;

    for (let i = 0n; i < s - 1n; i++) {
        x = bigIntPower(x, 2n, n);
        if (x === n - 1n) return true;
    }

    return false; // Composite
}

// Probabilistic primality test using Miller-Rabin
// Returns true if n is likely prime, false if certainly composite.
// k is the number of bases to test.
// For practical purposes in factorization context, 15-20 iterations are usually enough
// for numbers up to ~2^256 or more, yielding a very low probability of error.
function isPrime(n, k = 20) {
    if (n <= 1n) return false;
    if (n <= 3n) return true;
    if (n % 2n === 0n || n % 3n === 0n) return false;

    // For very small n, specific deterministic bases exist.
    // Given where this is called (after limited TD on potential large composite),
    // n is likely large. This check is a safeguard for smaller remainders.
     if (n < 50n) { // A small threshold for potentially deterministic test
         const smallDeterministicBases = [2n, 3n, 5n, 7n, 11n, 13n]; // Sufficient for N < 2*10^10 approx
         for (let a of smallDeterministicBases) {
             if (n === a) return true; // n is one of the bases
             if (n > a && !millerRabinTest(n, a)) return false;
         }
         return true; // Likely prime for this small range and bases
     }


    // For larger n, use k random bases in [2, n-2].
    const maxBase = n - 2n;
    const minBase = 2n;

    // If n is very small, n-2 < 2, cannot pick bases in [2, n-2].
    // This case is handled by n <= 3n check or the n < 50n check above.
    if (minBase > maxBase) {
        // This happens if n < 4. The checks above handle n<=3.
        // If n=4, minBase=2, maxBase=2. range is [2,2]. Handled by n%2==0.
        // This conditional should ideally not be reached for n > 3.
        // As a fallback, if n > 3 but somehow minBase > maxBase (e.g., n=4),
        // and it wasn't caught by previous checks, it implies n must be 4, which is composite.
        return false;
    }

    for (let i = 0; i < k; i++) {
        // Generate random base 'a' in [2, n-2]
        let a = randomBigInt(minBase, maxBase);
        if (!millerRabinTest(n, a)) {
            return false; // Certainly composite
        }
    }

    return true; // Likely prime
}

// Pollard's Rho (Brent's method) with Batched GCD optimization for finding *one* factor
// Returns a factor > 1 or 1n if no factor is found within the iteration limit or attempts,
// which could indicate the number is prime or hard for this method/limit.
// This function is now called *after* a primality test for efficiency on large numbers.
function pollardsRho(n, max_iterations = 200000n) { // Increased default iterations slightly
    // Basic checks - mostly redundant if called after isPrime and initial TD, but kept for robustness.
    // These checks are often done *before* calling Rho in a composite testing phase.
    if (n <= 1n) return 1n;
    if (n === 2n || n === 3n) return n;
    if (n % 2n === 0n) return 2n;
    if (n % 3n === 0n) return 3n;

    // Pollard's Rho attempts with different start points and constants
    // These constants can impact performance for specific numbers.
    const attempts = [
        { start_x: 2n, c: 1n },
        { start_x: 3n, c: 2n },
        { start_x: 5n, c: 3n }
    ];

    // Batch size for accumulating differences before GCD calculation.
    const BATCH_SIZE = 40n;

    for (let attempt of attempts) {
        let x = attempt.start_x;
        let c = attempt.c;

        // Ensure start_x and c are less than n and c is not 0.
        // For large n, x % n and c % n are typically x and c.
        if (n > 0n) {
            x = x % n;
            c = c % n;
        }
         // If n is small, x or c might become 0 or 1.
         // If x becomes 0 or 1, f(x) becomes predictable.
         // If c becomes 0, f(x) = x^2 % n, which might find factors related to squares.
         // If x==c==0, f(x)=0. If x=0, c!=0, f(0)=c. If x!=0, c=0, f(x)=x^2.
         // Using fixed small starts (2,3,5) is usually fine for large n.
        if (c === 0n) c = 1n; // Avoid f(x)=x^2 % n unless c=1 is also chosen

        let f = (val) => {
             let res = val * val;
             res += c;
             res %= n;
             return res;
        };

        let y = x;
        let l = 1n; // Segment length
        let steps_taken = 0n;

        let prod = 1n; // Batched product of differences
        let steps_since_last_gcd = 0n;

        // Brent's cycle detection with batched GCD
        while (steps_taken < max_iterations) {
            y = x; // Mark segment start

            // Advance x (hare) l steps, batching differences
            for (let i = 0n; i < l && steps_taken < max_iterations; ++i) {
                x = f(x); // Hare takes one step
                steps_taken++;

                let diff = x - y;
                // Ensure diff is positive modulo n
                diff = (diff % n + n) % n;

                if (diff === 0n) {
                    // Found x === y (mod n). This implies a cycle mod n, often means attempt failed.
                    steps_taken = max_iterations; // Signal failure for this attempt
                    break;
                }

                prod = (prod * diff) % n;

                steps_since_last_gcd++;

                // Perform batched GCD check periodically or at the end of segment/iterations
                if (steps_since_last_gcd === BATCH_SIZE || i === l - 1n || steps_taken >= max_iterations) {
                     // Only check if product is potentially useful (not 0 or 1 mod n).
                     if (prod > 1n) {
                        let d = bigIntGCD(n, prod);

                        if (d !== 1n) {
                            if (d === n) {
                                // Batched product is a multiple of n -> cycle mod n detected within batch
                                steps_taken = max_iterations; // Signal failure
                                break; // Break inner loop
                            }
                            // Found a non-trivial factor d.
                            return d; // Return the factor immediately
                        }
                     }
                     // If d was 1n (or prod was 0 or 1), reset batch.
                     // If d was n, we broke the inner loop.
                     if (steps_taken < max_iterations) { // Only reset if not already signaled failure
                         prod = 1n;
                         steps_since_last_gcd = 0n;
                     }
                }
                 if (steps_taken >= max_iterations) break; // Exit inner loop if max iterations reached
            }

            // After advancing l steps (or hitting max_iterations or finding a factor)
            // If a factor was found and returned, this code isn't reached.
            // If max_iterations was reached or a cycle mod n was detected (steps_taken = max_iterations), break outer.
            if (steps_taken >= max_iterations) break;

            // Double the segment length for the next Brent iteration
            l *= 2n;
        }
    }

    return 1n; // Indicate failure to find factor within limits
}

// Threshold for switching from trial division to Pollard's Rho/Miller-Rabin
// If the potential trial divisor exceeds this value, trial division stops
// and Rho/MR is attempted on the remaining number.
const TRIAL_DIVISION_LIMIT = 1_000_000n; // Factors up to 1M are found by trial division

/**
 * Finds the prime factors of a given positive integer (BigInt).
 * This version combines optimized trial division for small factors with
 * Miller-Rabin primality test and Pollard's Rho for larger numbers.
 * It's optimized for large numbers by using MR to avoid Rho on primes.
 *
 * @param {bigint} number The positive integer (as BigInt) to factorize.
 * @returns {bigint[]} An array containing the prime factors of the number,
 *                     including duplicates if a factor appears multiple times.
 *                     Returns an empty array for numbers less than 2.
 */
function primeFactorization(number) {
    // Ensure input is BigInt
    if (typeof number !== 'bigint') {
        try {
            number = BigInt(number);
        } catch (e) {
            throw new Error("Input must be a BigInt or coercible to BigInt.");
        }
    }

    if (number < 2n) {
        return [];
    }
    // Basic base cases for small primes already handled by trial division, but explicit is clear.
    // if (number === 2n) return [2n]; // Handled by initial TD
    // if (number === 3n) return [3n]; // Handled by initial TD

    const factors = [];
    let n = number; // Use n internally

    // --- Initial Trial Division for 2 and 3 ---
    while (n % 2n === 0n) {
        factors.push(2n);
        n /= 2n;
    }
    while (n % 3n === 0n) {
        factors.push(3n);
        n /= 3n;
     }
     if (n === 1n) {
        // If only 2s and 3s were factors, sorting is optional but good practice.
        // The main sort at the end handles this too, but sorting early might be slightly faster if the stack is empty.
        factors.sort((a, b) => a < b ? -1 : (a > b ? 1 : 0));
        return factors;
     }

    // Use a stack to manage numbers that need factoring after initial TD by 2/3.
    const numbersToFactorStack = [n];
    let numToFactor; // Variable to hold the current number being factored from the stack

    // Pre-calculate the square of the trial division limit.
    // Numbers <= this limit squared can be fully factored by trial division alone.
    const TD_LIMIT_SQUARED = TRIAL_DIVISION_LIMIT * TRIAL_DIVISION_LIMIT;

    // Number of bases to test for Miller-Rabin.
    const MILLER_RABIN_ITERATIONS = 20;

    // Process numbers from the stack until it's empty
    while ((numToFactor = numbersToFactorStack.pop()) !== undefined) {
        if (numToFactor === 1n) continue; // Already factored

        // Optimization for smaller numbers: If the number is small enough
        // (numToFactor <= TD_LIMIT_SQUARED), use full trial division up to its square root.
        // This also correctly identifies primes within this range efficiently.
        if (numToFactor <= TD_LIMIT_SQUARED) {
             // Apply standard trial division up to sqrt(numToFactor)
             let trial_i = 5n;
             let limit_sqrt_num = numToFactor; // Use numToFactor itself for the sqrt check limit

             // Loop while the square of the current trial divisor is <= the number being factored.
             while (trial_i * trial_i <= limit_sqrt_num) {
                 // Check factor i (6k - 1)
                 while (numToFactor % trial_i === 0n) {
                     factors.push(trial_i);
                     numToFactor /= trial_i;
                 }
                 if (numToFactor === 1n) break; // Fully factored this number

                 let trial_i_plus_2 = trial_i + 2n;

                 // Check factor i+2 (6k + 1)
                 while (numToFactor % trial_i_plus_2 === 0n) {
                     factors.push(trial_i_plus_2);
                     numToFactor /= trial_i_plus_2;
                 }
                 if (numToFactor === 1n) break; // Fully factored this number

                 trial_i += 6n; // Move to next pair (i+6)
             }
             // After full trial division up to sqrt, if numToFactor > 1, the remainder is prime.
             if (numToFactor > 1n) {
                factors.push(numToFactor);
             }
             continue; // Finished factoring this number from stack
        }

        // --- Hybrid Factoring for numToFactor (when it's large, i.e., > TD_LIMIT_SQUARED) ---
        // For numbers larger than TD_LIMIT_SQUARED, their smallest prime factor might be
        // larger than TD_LIMIT. Apply trial division only up to TD_LIMIT.
        let trial_i = 5n;
        const TD_LIMIT = TRIAL_DIVISION_LIMIT; // Use constant for clarity

        // Trial division up to TRIAL_DIVISION_LIMIT.
        // We already handled cases where sqrt <= TD_LIMIT above.
        // So here, sqrt(numToFactor) > TD_LIMIT. We only need to check trial divisors up to TD_LIMIT.
        while (trial_i <= TD_LIMIT) { // Loop up to the limit itself
             // Check factor i (6k - 1)
            while (numToFactor % trial_i === 0n) {
                factors.push(trial_i);
                numToFactor /= trial_i;
            }
            if (numToFactor === 1n) break; // Fully factored this number

            let trial_i_plus_2 = trial_i + 2n;

            // Check factor i+2 (6k + 1)
            // Only check if i+2 is also within the trial division limit.
            if (trial_i_plus_2 > TD_LIMIT) {
                 // The next candidate (i+2) exceeds the limit. Stop trial division *after* checking trial_i.
                break;
            }

            while (numToFactor % trial_i_plus_2 === 0n) {
                factors.push(trial_i_plus_2);
                numToFactor /= trial_i_plus_2;
            }
            if (numToFactor === 1n) break; // Fully factored this number

            trial_i += 6n; // Move to next pair (i+6)
        }
        // End of Limited Trial Division Loop

        // After limited trial division, if numToFactor > 1, its smallest prime factor is > TD_LIMIT.
        if (numToFactor === 1n) {
            continue; // This number was fully factored by limited trial division
        }

        // --- Primality Test BEFORE Pollard's Rho ---
        // If the remaining numToFactor is likely prime, add it directly.
        // This is the core optimization for large numbers.
        if (isPrime(numToFactor, MILLER_RABIN_ITERATIONS)) {
            // numToFactor is likely prime. Add it to factors.
            factors.push(numToFactor);
        } else {
            // numToFactor is composite. Attempt Pollard's Rho to find a factor.
            let factor = pollardsRho(numToFactor);

            if (factor > 1n && factor < numToFactor) {
                // Rho found a non-trivial factor. Add the factor and the remaining quotient to the stack
                // to be factored further.
                numbersToFactorStack.push(numToFactor / factor);
                numbersToFactorStack.push(factor);
            } else {
                // Rho failed (returned 1n). This happens if Rho couldn't find a factor
                // within the iteration limits or given the starting parameters.
                // Since we already know it's composite from Miller-Rabin, adding
                // the remaining numToFactor as a prime factor here is logically incorrect.
                // A robust factorizer would try other algorithms (like ECM, MPQS) here
                // or report failure to fully factor.
                // However, following the original algorithm's structure on Rho failure
                // (which was to add the unfactored number as if it were prime), we keep this.
                // This indicates a limitation that this number couldn't be factored
                // by the combined TD+Rho methods within limits, despite being composite.
                 factors.push(numToFactor); // Fallback: add unfactored composite
            }
        }
    } // End of stack processing loop

    // Sort factors for canonical representation (optional but good for testing)
    factors.sort((a, b) => a < b ? -1 : (a > b ? 1 : 0));

    return factors;
}

// Main function to test the primeFactorization function
function main() {
    const testCases = [
        0n,
        1n,
        2n,
        3n,
        4n,
        6n,
        25n,
        30n,
        72n, // 2*2*2 * 3*3
        97n, // A prime number
        100n,
        13195n, // Project Euler problem 3 example: 5 * 7 * 13 * 29
        600851475143n, // Project Euler problem 3 number
        2n ** 32n, // A large power of 2
        3n ** 20n, // A large power of 3
        (2n ** 10n) * (3n ** 5n) * (5n ** 3n) * 7n, // Composite with several small primes
        1000000007n, // A moderately large prime (10^9 + 7)
        1000000007n * 2n, // Product with 2
        1000000007n * 3n, // Product with 3
        1000000007n * 5n, // Product with 5
        1000000007n * 1000000007n, // Square of a prime
        1000000007n * 1000000009n, // Product of two large primes (10^9 + 7 and 10^9 + 9; note 10^9+9 is prime)
                                    // 1000000009 is actually prime.
        999999999989n, // A large prime (Largest prime < 10^12)
        (2n ** 60n), // Very large power of 2 to test early exit
        12345678901234567890n, // A very large composite number
        // A number designed to test the (i+2)*(i+2) > number break condition more directly
        // if the remaining number is prime. E.g., a prime just above a square.
        // 47 is prime. i=5. 5*5=25 <= 47. i+2=7. 7*7=49 > 47. Break. Add 47.
        47n,
        // Product of two primes p1, p2 where p1 is 6k-1 and p2 is 6k+1, e.g., 29*31
        29n * 31n, // 899
        // A number that's just a large prime, to test the final 'if (number > 1n)'
        BigInt("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000067"), // 10^100 + 67 (a known large prime)
        // A number that is a product of a very small prime and a very large prime
        3n * BigInt("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000067"),
        
        BigInt("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000067"),
        
        10000000000000000000000000000000000000039n,
        
        100000000000000000000000005280000000000000000000000000066671n,
        
        100000000000000000000003n * 100000000000000000000005n,
        
        111111111111111111111111113n * 111111111111111111111111123n,
        
        100000000000000000000000000104000000000000000000000000002623n,
        
        1101996010049018383040097869069000000000000000000000000013777n
    ];

    console.log("Starting prime factorization tests...\n");

    for (const num of testCases) {
        console.log(`Input: ${num}`);
        const startTime = process.hrtime.bigint(); // For more precise timing
        let factors;
        let errorOccurred = false;
        try {
            factors = primeFactorization(num);
            const endTime = process.hrtime.bigint();
            const durationMs = Number(endTime - startTime) / 1e6; // Convert nanoseconds to milliseconds

            console.log(`Factors: [${factors.join(', ')}]`);
            console.log(`Time taken: ${durationMs.toFixed(3)} ms`);

            // Verification
            if (num < 2n) {
                if (factors.length === 0) {
                    console.log("Verification: PASSED (Correctly empty for num < 2n)");
                } else {
                    console.error("Verification: FAILED (Should be empty for num < 2n)");
                }
            } else if (factors.length === 0 && num >=2n) {
                 console.error(`Verification: FAILED (Resulted in empty array for ${num} >= 2n)`);
            }
            else {
                const product = factors.reduce((acc, val) => acc * val, 1n);
                if (product === num) {
                    console.log("Verification: PASSED (Product of factors matches original number)");
                } else {
                    console.error(`Verification: FAILED (Product ${product} !== Original ${num})`);
                }
            }
        } catch (e) {
            errorOccurred = true;
            const endTime = process.hrtime.bigint();
            const durationMs = Number(endTime - startTime) / 1e6;
            console.error(`Error for ${num}: ${e.message}`);
            console.log(`Time taken before error: ${durationMs.toFixed(3)} ms`);
        }
        console.log("------------------------------------");
    }
    console.log("\nPrime factorization tests finished.");
}

// Run the main function
// This check ensures main() runs only when the script is executed directly
if (require.main === module || (typeof process !== 'undefined' && process.env.NODE_ENV === 'test')) {
    main();
}
