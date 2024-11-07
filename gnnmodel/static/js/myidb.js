const dbname = "myidatabase",
  storename = "gnnepcsaftparameters";

async function makingdb() {
  const db = await idb.openDB(dbname, 1, {
    upgrade(db) {
      const gnnstore = db.createObjectStore(storename, {
        keyPath: "id",
        autoIncrement: true,
      });

      gnnstore.createIndex("inchi", "inchi", { unique: true });
    },
  });
  alasql(
    'ATTACH INDEXEDDB DATABASE myidatabase; \
      USE myidatabase; \
      SELECT * INTO gnnepcsaftparameters FROM CSV("/static/mydata.csv", {headers:true, separator: " | "});',
    [],
    function (res) {
      console.log("Number of records added: ", res[2]);
    }
  );
  return db;
}

const db = makingdb();

async function getinchi(key) {
  return (await db).getFromIndex(storename, "inchi", key);
}
