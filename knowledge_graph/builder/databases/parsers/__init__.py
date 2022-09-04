from builder.databases.parsers.cancer_genome_interpreter_parser import CancerGenomeInterpreterParser
from builder.databases.parsers.drug_bank_parser import DrugBankParser
from builder.databases.parsers.uniprot_parser import UniProtParser
from builder.databases.parsers.hgnc_parser import HGNCParser
from builder.databases.parsers.hmdb_parser import HMDBParser
from builder.databases.parsers.reactome_parser import ReactomeParser
from builder.databases.parsers.foodb_parser import FooDBParser


parsers = {
    "CancerGenomeInterpreter": CancerGenomeInterpreterParser,
    "DrugBank": DrugBankParser,
    "UniProt": UniProtParser,
    "HGNC": HGNCParser,
    "HMDB": HMDBParser,
    "Reactome": ReactomeParser,
    "FooDB": FooDBParser
}
