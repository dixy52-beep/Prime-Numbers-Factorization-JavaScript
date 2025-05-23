// Helper function for Greatest Common Divisor (GCD) using Euclidean algorithm for BigInt
function bigIntGCD(a, b) {
    // Handle negative inputs by taking absolute value.
    a = a < 0n ? -a : a;
    b = b < 0n ? -b : b;
    while (b) {
        [a, b] = [b, a % b];
    }
    return a;
}

// Pollard's Rho (Brent's method) with Batched GCD optimization for finding *one* factor
// Returns a factor > 1 or 1n if no factor is found within the iteration limit or attempts,
// which could indicate the number is prime or hard for this method/limit.
// Note: A robust primality test (like Miller-Rabin) would be needed for certainty,
// but is skipped here for simplicity in demonstrating the algorithm switch and Rho optimization.
function pollardsRho(n, max_iterations = 200000n) { // Increased default iterations slightly
    if (n <= 1n) return 1n; // Should not happen if called correctly
    if (n === 2n || n === 3n) return n; // Return the number itself if it's 2 or 3 (already prime)
    if (n % 2n === 0n) return 2n; // Should be handled by trial division, but as a safeguard
    if (n % 3n === 0n) return 3n; // As a safeguard

    // The small prime check is mostly handled by the main factorization function's trial division.

    // Pollard's Rho attempts with different start points and constants
    // Using a few fixed attempts. More robust versions might use random choices.
    const attempts = [
        { start_x: 2n, c: 1n },
        { start_x: 3n, c: 2n },
        { start_x: 5n, c: 3n } // Changed one starting point for variety
    ];

    // Batch size for accumulating differences before GCD calculation.
    // Larger size reduces GCD calls but increases product computation cost.
    // A value around 20-50 is typical. Choosing 40n.
    const BATCH_SIZE = 40n;

    for (let attempt of attempts) {
        let x = attempt.start_x;
        let c = attempt.c;

        // Ensure start_x and c are less than n and c is not 0 (usually).
        // If n is very small, taking modulo might make x or c zero, which is handled below.
        if (n > 0n) { // Avoid modulo by zero or negative
            x = x % n;
            c = c % n;
        }
        // If after modulo c becomes 0, change it to 1 to avoid trivial f(x)=x*x.
        if (c === 0n) c = 1n;
        // If after modulo x becomes 0, 1, 2, this might be less effective.
        // For simplicity here, we use the fixed start points. For robustness,
        // ensure start_x is not a trivial value like 0, 1, n-1, etc.
        // The initial fixed points (2,3,5) are generally fine unless n is very small.

        // The pseudo-random sequence function: f(val) = (val*val + c) % n
        let f = (val) => {
             let res = val * val;
             res += c;
             res %= n;
             return res;
        };

        let y = x; // Tortoise (slow pointer)
        let l = 1n; // Current segment length (powers of 2)
        let steps_taken = 0n; // Total steps taken for this attempt

        let prod = 1n; // Product of differences (x_k - y) % n within the current batch
        let steps_since_last_gcd = 0n;

        // Brent's cycle detection with batched GCD
        // The outer loop structure advances x by l steps, checks GCD, doubles l.
        while (steps_taken < max_iterations) {
            y = x; // y is the start of the current segment (where x was before advancing l steps)

            // Advance x (hare) l steps
            // The inner loop advances x step by step, accumulating product for batching.
            for (let i = 0n; i < l && steps_taken < max_iterations; ++i) {
                x = f(x); // Hare takes one step
                steps_taken++;

                let diff = x - y; // Difference relative to segment start (y)

                // Ensure diff is positive modulo n before multiplying.
                // Using (diff % n + n) % n handles both positive and negative results of diff % n
                // correctly in JavaScript for BigInts.
                diff = (diff % n + n) % n;

                // If diff is 0, we found a cycle mod n instantly (x == y mod n).
                // This implies the attempt failed.
                if (diff === 0n) {
                    steps_taken = max_iterations; // Signal failure for this attempt
                    break; // Exit inner loop
                }

                prod = (prod * diff) % n; // Accumulate product modulo n

                steps_since_last_gcd++;

                // Perform batched GCD check periodically or at the end of the l steps
                // Check when a batch is complete OR if we are at the last step of the current l-segment (i == l-1n)
                // Also check if we are about to exceed max_iterations after this step.
                if (steps_since_last_gcd === BATCH_SIZE || i === l - 1n || steps_taken >= max_iterations) {
                     // Only check if product is potentially useful (not 0 or 1 mod n).
                     if (prod > 1n) {
                        let d = bigIntGCD(n, prod);

                        if (d !== 1n) {
                            // Found a factor related to this batch's product.
                            // If d < n, it's a non-trivial factor. Return it.
                            // If d === n, the batched product is a multiple of n. This suggests a cycle mod n.
                            if (d === n) {
                                // Batched product is a multiple of n. This means one of the (x_k - y) terms was a multiple of n,
                                // or the product of several (x_k - y) terms is a multiple of n.
                                // This typically indicates a cycle mod n was found within this batch.
                                // Treat as failure for this attempt.
                                steps_taken = max_iterations; // Signal failure for this attempt
                                break; // Break inner for loop
                            }
                            // Found a non-trivial factor d. This factor might be composite.
                            // The factorization function handles factoring composite factors.
                            return d; // Return the factor immediately
                        }
                     }
                     // If d was 1n (or prod was 0 or 1), reset batch product and counter
                     // If d was n, we broke.
                     if (steps_taken < max_iterations) { // Only reset if not already signaled failure
                         prod = 1n;
                         steps_since_last_gcd = 0n;
                     }
                } // End of batch check

                 if (steps_taken >= max_iterations) break; // Exit inner loop if max iterations reached
            } // End of for loop (l steps advance)

            // After advancing l steps (or hitting max_iterations or finding a factor)
            // If a factor was found and returned, this code isn't reached.
            // If max_iterations was reached or a cycle mod n was detected (steps_taken = max_iterations), break outer.
            if (steps_taken >= max_iterations) break;

            // Double the segment length for the next Brent iteration
            l *= 2n;

            // Add a safeguard for excessively large l relative to steps_taken?
            // max_iterations should handle overall progress limit.
        } // End of while loop (Brent's)

        // If while loop exited without returning a factor (max_iterations reached or cycle mod n)
        // No factor found for this attempt within the limits.
    } // End of attempts loop

    // If all attempts failed to find a factor within the limit
    return 1n; // Indicate failure
}

// Threshold for switching from trial division to Pollard's Rho
// If the potential trial divisor exceeds this value, trial division stops
// and Pollard's Rho is attempted on the remaining number.
// Set this based on what's feasible for trial division on typical hardware.
// 1,000,000 is a common balance point. sqrt(10^12) = 1M.
// This means trial division is used for factors up to 1M, effectively handling
// numbers up to ~10^12 quite quickly. For larger numbers, Rho is used.
const TRIAL_DIVISION_LIMIT = 1_000_000n;

/**
 * Finds the prime factors of a given positive integer (BigInt).
 * This version combines optimized trial division for small factors with
 * Pollard's Rho for larger composite numbers. It is designed to be
 * significantly faster than pure trial division for numbers with large
 * prime factors that are not excessively large (e.g., up to ~40-50 digits).
 *
 * This version includes batched GCD checks within Pollard's Rho for improved
 * performance during the Rho phase and refines the trial division application
 * based on the size of the number being factored.
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

    const factors = [];
    let n = number; // Use n internally

    // --- Trial Division for 2 ---
    while (n % 2n === 0n) {
        factors.push(2n);
        n /= 2n;
    }
    if (n === 1n) return factors;

    // --- Trial Division for 3 ---
    while (n % 3n === 0n) {
        factors.push(3n);
        n /= 3n;
     }
     if (n === 1n) return factors;

    // Use a stack to manage numbers that need factoring.
    // Start with the remaining number after handling 2 and 3.
    const numbersToFactorStack = [n];
    let numToFactor; // Variable to hold the current number being factored from the stack

    // Pre-calculate the square of the trial division limit.
    // Numbers <= this limit squared can be fully factored by trial division alone.
    const TD_LIMIT_SQUARED = TRIAL_DIVISION_LIMIT * TRIAL_DIVISION_LIMIT;

    // Process numbers from the stack until it's empty
    while ((numToFactor = numbersToFactorStack.pop()) !== undefined) {
        if (numToFactor === 1n) continue; // Already factored

        // Optimization: If the number on the stack is small enough that its sqrt
        // is less than or equal to the TRIAL_DIVISION_LIMIT, just use trial division
        // fully up to its square root on it. This avoids calling Rho unnecessarily
        // on smaller composite numbers where trial division is already efficient.
        // Condition: sqrt(numToFactor) <= TD_LIMIT  <=>  numToFactor <= TD_LIMIT * TD_LIMIT
        if (numToFactor <= TD_LIMIT_SQUARED) {
             // Apply standard trial division up to sqrt(numToFactor)
             let trial_i = 5n;
             let limit_sqrt_num = numToFactor; // Use numToFactor itself for the sqrt check limit

             // Loop while the square of the current trial divisor is <= the number being factored.
             // This is the standard O(sqrt(N)) trial division limit.
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
             continue; // Finished factoring this number from the stack
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

        // After limited trial division:
        // - If numToFactor is 1, it was fully factored by trial division up to the limit.
        // - If numToFactor > 1, its smallest prime factor is > TD_LIMIT.
        //   In this case, it must be factored using Pollard's Rho or similar methods.

        if (numToFactor === 1n) {
            continue; // This number is fully factored
        }

        // If trial division didn't fully factor the number (numToFactor > 1n),
        // attempt Pollard's Rho on the remaining numToFactor.
        // Note: If numToFactor is actually prime and > TD_LIMIT_SQUARED, pollardsRho might return 1n.
        // A primality test would be ideal here before calling Rho, but skipped per original structure.
        let factor = pollardsRho(numToFactor);

        if (factor > 1n && factor < numToFactor) {
            // Rho found a non-trivial factor. Add the factor and the remaining quotient to the stack
            // to be factored further.
            // Push the larger number first to potentially process a larger chunk sooner,
            // though order doesn't strictly matter for correctness with a stack/queue.
            numbersToFactorStack.push(numToFactor / factor);
            numbersToFactorStack.push(factor);
        } else {
            // Rho failed (returned 1n) or found a trivial factor (shouldn't happen for n > 3n).
            // The remaining numToFactor is likely prime (especially if > TD_LIMIT_SQUARED)
            // or too hard for this Rho implementation within the limit.
            // In a real system, you'd try ECM or other methods here.
            // For this simplified version, add the remaining numToFactor as a prime factor.
            // This acts as the fallback for prime remainders or Rho failures.
            factors.push(numToFactor);
        }
    }

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
        3n * BigInt("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000067")
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
