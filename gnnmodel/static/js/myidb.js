const dbname = "gnnepcsaft_idb",
  storename = "gnnepcsaftparameters";

async function makingdb() {
  const db = await idb.openDB(dbname, 1, {
    upgrade(db) {
      const gnnstore = db.createObjectStore(storename, {
        keyPath: "inchi",
      });

      gnnstore.createIndex("smiles", "smiles", { unique: false });

      alasql(
        `ATTACH INDEXEDDB DATABASE ${dbname}; 
        USE ${dbname}; 
        SELECT ROUND(m, 4) as m, ROUND(sigma, 4) as sigma, ROUND(e, 4) as e, 
        ROUND(k_ab, 4) as k_ab, ROUND(e_ab, 4) as e_ab, mu, na, nb, inchi, smiles 
        INTO ${storename} 
        FROM CSV("/static/mydata.csv", {headers:true, separator: "|"});`,
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
