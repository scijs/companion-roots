var findRoots = require("../roots.js")

require("tap").test("companion roots", function(t) {

  //Real roots
  var rroots = findRoots([1, 1, -1])
  
  //Complex roots
  var croots = findRoots([1, 1, -1], [0, 0, 0])
  
  t.end()
})
