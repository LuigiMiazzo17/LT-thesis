# Ottimizzazione della libreria CoopeRIS tramite calcolo parallelo

[![ci](https://github.com/unitn-drive/thesis/actions/workflows/ci.yaml/badge.svg)](https://github.com/unitn-drive/thesis/actions/workflows/ci.yaml)

## Abstract

Le Reconfigurable Intelligent Surfaces (RISs) sono una tecnologia emergente che \
promette di rivoluzionare le comunicazioni wireless del futuro, permettendo di \
modificare in tempo reale le caratteristiche del canale di comunicazione tra un \
trasmettitore e un ricevitore quando non è possibile instaurare una connessione \
line-of-sight. Tuttavia, lo studio e la valutazione della fattibilità e delle \
prestazioni di questa tecnologia, come per molte altre, richiede l'utilizzo di \
strumenti di simulazione sofisticati e sviluppati ad hoc. I framework di \
simulazione sono essenziali per lo studio di tecnologie ancora in uno stato \
embrionale come le RIS, poiché permettono di testare e confrontare diverse \
configurazioni e algoritmi in un ambiente controllato e riproducibile, evitando \
costosi test sperimentali, in particolar modo per prodotti non ancora \
disponibili sul mercato. Talvolta però, la complessità dei modelli simulativi e \
la mole di dati da processare richiedono un'elevata potenza di calcolo e \
risorse hardware, che possono risultare insufficienti per ottenere risultati in \
tempi ragionevoli. In questo contesto, il calcolo parallelo offre la \
possibilità di sfruttare al meglio le risorse a disposizione, riducendo i tempi \
di esecuzione e migliorando le prestazioni, tuttavia complicando in alcuni casi \
la portabilità e la manutenibilità del codice. Questo elaborato si propone di \
analizzare tre diversi framework di programmazione parallela confrontando le \
principali caratteristiche che li contraddistinguono, nell'ottica di valutarne \
l'efficacia, la versatilità, e ove applicabile la difficoltà d'integrazione. I \
framework presi in considerazione sono: parallelizzazione su CPU tramite la \
libreria libpthread e su GPU tramite le librerie CUDA e OpenCL. Per questo \
scopo, è stata eseguita l'integrazione di ciascuno di questi paradigmi nella \
libreria CoopeRIS, framework di simulazione per le RISs, al fine di valutarne \
l'impatto sulle prestazioni delle simulazioni. I risultati ottenuti da questo \
studio dimostrano pienamente l'efficacia del calcolo parallelo, in particolare \
su GPU, nel ridurre i tempi di esecuzione delle citate simulazioni, \
essenzialmente indicando un miglioramento pressoché lineare nella velocità di \
esecuzione al crescere del numero di thread utilizzati in riferimento \
all'implementazione su CPU, e un incremento di quasi due ordini di grandezza \
tramite le implementazioni su GPU.

## License

This project is licensed under the [MIT](https://opensource.org/licenses/MIT) \
license See [LICENSE](./LICENSE) file for details
