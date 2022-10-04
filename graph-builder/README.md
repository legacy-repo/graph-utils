# Knowledge Graph
A tool for building a knowledge graph by parsing and connecting more than twenty databases.

```bash
(biomedgps) ➜ /Users/codespace/Documents/Code/BioMedGPS git:(master) ✗ > graph-builder --help
Usage: graph-builder [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  parse-database  Parse databases and make the related graph files.
  parse-ontology  Parse ontologies and make related graph files.
  print-config   Print the default config file.
```

### How to use it?

#### Step 1. Download and parse ontology databases

```bash
(biomedgps) ➜ /Users/codespace/Documents/Code/BioMedGPS git:(master) ✗ > graph-builder parse-ontology --help
Usage: graph-builder parse-ontology [OPTIONS]

  Parse ontologies and make related graph files.

Options:
  -d, --ontology-dir PATH     The directory which saved the downloaded
                              database files.  [required]

  -o, --output-dir PATH       The directory which saved the graph files.
                              [required]

  --download / --no-download  Whether download the source file(s)?
  --skip / --no-skip          Whether skip the existing file(s)?
  --help                      Show this message and exit.
```

#### Step 2. Download and parse databases

```bash
(biomedgps) ➜ /Users/codespace/Documents/Code/BioMedGPS git:(master) ✗ > graph-builder parse-database --help
Usage: graph-builder parse-database [OPTIONS]

  Parse databases and make the related graph files.

Options:
  -d, --db-dir PATH               The directory which saved the downloaded
                                  database files.  [required]

  -o, --output-dir PATH           The directory which saved the graph files.
                                  [required]

  -c, --config PATH               The config file related with database.
  --database [DrugBank|UniProt|HGNC|HMDB|Reactome|FooDB|CancerGenomeInterpreter|PathwayCommons|OncoKB|Pfam|RefSeq|DisGEnet|STITCH|STRING|IntAct|MutationDs|JensenLab|GWASCatalog|DGIdb|SIDER|PhosphoSitePlus|SIGNOR|CORUM|ExposomeExplorer|HPA|SMPDB]
                                  Which databases (you can specify the
                                  --database argument multiple times)?
                                  [required]

  -n, --n-jobs INTEGER            Hom many jobs?
  --download / --no-download      Whether download the source file(s)?
  --skip / --no-skip              Whether skip the existing file(s)?
  --help                          Show this message and exit.
```