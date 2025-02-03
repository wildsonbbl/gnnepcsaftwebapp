const dbname = "myidatabase",
  storename = "gnnepcsaftparameters";

async function makingdb() {
  const db = await idb.openDB(dbname, 1, {
    upgrade(db) {
      const gnnstore = db.createObjectStore(storename, {
        keyPath: "inchi",
      });

      gnnstore.createIndex("smiles", "smiles", { unique: false });

      alasql(
        'ATTACH INDEXEDDB DATABASE myidatabase; \
         USE myidatabase; \
         SELECT m, sigma, e, k_a, e_ab, mu, na, nb, inchi, smiles \
         INTO gnnepcsaftparameters \
         FROM CSV("/static/mydata.csv", {headers:true, separator: "|"});',
        [],
        function (res) {
          console.log("Number of records added: ", res[2]);
        }
      );
    },
  });

  return db;
}

const db = makingdb();

async function getinchi(key) {
  return (await db).get(storename, key);
}

async function getsmiles(key) {
  return (await db).getFromIndex(storename, "smiles", key);
}
