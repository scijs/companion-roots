"use strict"

var numeric = require("numeric")
var max = Math.max
var abs = Math.abs
var EPSILON = 1e-10

function findRoots_cplx(r_coeff, i_coeff) {
  //Trim zero coefficients
  var N = r_coeff.length
    , r_c, i_c, w
  while(N-- > 0) {
    r_c = r_coeff[N]
    i_c = i_coeff[N]
    w = r_c*r_c + i_c*i_c
    if(w > EPSILON) {
      break
    }
  }
  
  //Handle degenerate case
  if(N === 0) {
    return [[], []]
  }
  
  //Assemble matrix
  var A = numeric.rep([2*N, 2*N], 0.0)
  for(var i=2; i<2*N; ++i) {
    A[i-2][i] = 1.0
  }
  
  //Compute negative reciprocal of leading coeff
  w = 1.0 / w
  r_c *= -w
  i_c *= w
  
  //Compute tensor product expansion of companion matrix
  var s_c = r_c + i_c
    , d_c = i_c - r_c
    , k1, k2, k3, rr, ii, r_a, i_a
    , R1 = A[2*N-2]
    , R2 = A[2*N-1]
  for(var i=0; i<N; ++i) {
    r_a = r_coeff[i]
    i_a = i_coeff[i]
    k1 = (r_a + i_a) * r_c
    k2 = r_a * d_c
    k3 = i_a * s_c
    rr = k1 - k3
    ii = k1 + k2
    R1[2*i] = R2[2*i+1] = rr
    R1[2*i+1] = -ii
    R2[2*i] = ii
  }
  
  //Compute roots
  var roots = numeric.eig(A).lambda
  var r_roots = roots.x
  var i_roots = roots.y || numeric.rep([2*N], 0.0)
  
  //Remove duplicate pairs of roots by brute force search
  var dupes = numeric.rep([2*N], false)
  var ptr = 0
  for(var i=0; i<2*N; ++i) {
    if(dupes[i]) {
      continue
    }
    var closest_idx = i+1
    var closest_d = Infinity
    for(var j=i+1; j<2*N; ++j) {
      if(dupes[j]) {
        continue
      }
      var d = max(abs(r_roots[i]-r_roots[j]), abs(i_roots[i]-i_roots[j]))
      if(d < closest_d) {
        closest_idx = j
        closest_d = d
      }
    }
    dupes[closest_idx] = true
    r_roots[ptr] = r_roots[i]
    i_roots[ptr] = i_roots[i]
    ++ptr
  }
  r_roots.length = N
  i_roots.length = N
  return [r_roots, i_roots]
}


function findRoots_real(coeff) {
  //Trim zero coefficients
  var N = coeff.length
  while(N > 0 && abs(coeff[--N]) < EPSILON) {
  }
  
  //Handle degenerate case
  if(N === 0) {
    return [[], []]
  }
  
  //Assemble matrix
  var A = numeric.rep([N, N], 0.0)
  for(var i=1; i<N; ++i) {
    A[i-1][i] = 1.0
  }
  var R = A[N-1]
  var c = -1.0 / coeff[N]
  for(var i=0; i<N; ++i) {
    R[i] = coeff[i] * c
  }
  
  //Compute roots
  var roots = numeric.eig(A).lambda
  if(roots.y) {
    return [ roots.x, roots.y ]
  } else {
    return [ roots.x, numeric.rep([N], 0.0) ]
  }
}

function findRoots(r_coeff, i_coeff) {
  if(i_coeff) {
    return findRoots_cplx(r_coeff, i_coeff)
  }
  return findRoots_real(r_coeff)
}

module.exports = findRoots