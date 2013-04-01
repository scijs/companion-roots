companion-roots
===============
Finds all [roots](http://en.wikipedia.org/wiki/Root_of_a_function) of a [polynomial](http://en.wikipedia.org/wiki/Polynomial) by computing the [eigenvalues](http://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors) of its [companion matrix](http://en.wikipedia.org/wiki/Companion_matrix).  In other words, it [factorizes the polynomial](http://en.wikipedia.org/wiki/Factorization_of_polynomials) over the complex numbers.

Use
===
Install using npm:

    npm install companion-roots
    
```javascript
var findRoots = require("companion-roots")

var roots = findRoots([1, 1, -1])  // Finds roots for 1 + 1*x - 1*x^2

// Now:
//      roots[0] = real part of roots
//      roots[1] = imaginary part of roots

for(var i=0; i<roots.length; ++i) {
  console.log(roots[0][i] + "+" + roots[1][i] + "i")
}

// Prints:
//  1.618033988749895+0i
//  -0.6180339887498949+0i
```

### `require("companion-roots")(real_coeffs[, imag_coeffs])`
Computes the roots of a polynomial

* `real_coeffs` the real coefficients of the polynomial arranged in order of increasing degree
* `imag_coeffs` (optional) the imaginary coefficients of the polynomial.  If not specified, assumed to be zero

**Returns:** A pair of vectors representing the real and imaginary parts of the roots of the polynomial

**Time Complexity:** `O(real_coeffs.length^3)`

Credits
=======
(c) 2013 Mikola Lysenko. MIT License