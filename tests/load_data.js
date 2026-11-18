// Load the generated data/*.js files in node (they call PW.data(kind, key, obj)).
const fs = require("fs"), path = require("path"), vm = require("vm");
const DATA = path.join(__dirname, "..", "docs", "data");
function load(file) {
  const store = {};
  const ctx = { PW: { data: (kind, key, obj) => { store[kind + ":" + key] = obj; } } };
  vm.runInNewContext(fs.readFileSync(path.join(DATA, file), "utf8"), ctx);
  return Object.values(store)[0];
}
module.exports = { load, DATA };
