import { openDB, deleteDB, wrap, unwrap } from 'https://unpkg.com/idb?module';

const dbname = "myidatabase",
      storename = "gnnepcsaftparameters";

const db = await openDB(dbname, 1, {
  upgrade(db) {
    
    const gnnstore = db.createObjectStore(storename, {
      keyPath: 'id',
      autoIncrement: true,
    });

    gnnstore.createIndex('inchi', 'inchi', {unique: true});
  },
});

// instância da transação
const store = db.transaction(storename, 'readwrite').objectStore(storename);



// operações assíncronas a serem resolvidas
await Promise.all([
 store.add({
  inchi: 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
  m: 5.87,
  sigma: 2.31,
  e:  187.87
 }),
]);