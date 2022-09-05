from builder.databases.parsers.cancer_genome_interpreter_parser import CancerGenomeInterpreterParser
from builder.databases.parsers.drug_bank_parser import DrugBankParser
from builder.databases.parsers.uniprot_parser import UniProtParser
from builder.databases.parsers.hgnc_parser import HGNCParser
from builder.databases.parsers.hmdb_parser import HMDBParser
from builder.databases.parsers.reactome_parser import ReactomeParser
from builder.databases.parsers.foodb_parser import FooDBParser
from builder.databases.parsers.pathway_commons_parser import PathwayCommonsParser
from builder.databases.parsers.oncokb_parser import OncoKBParser
from builder.databases.parsers.pfam_parser import PfamParser
from builder.databases.parsers.refseq_parser import RefSeqParser
from builder.databases.parsers.disgenet_parser import DisGEnetParser
from builder.databases.parsers.stitch_parser import STITCHParser
from builder.databases.parsers.string_parser import STRINGParser
from builder.databases.parsers.intact_parser import IntActParser
from builder.databases.parsers.mutationds_parser import MutationDsParser
from builder.databases.parsers.jensenlab_parser import JensenLabParser
from builder.databases.parsers.gwas_catalog_parser import GWASCatalogParser
from builder.databases.parsers.dgi_db_parser import DGIdbParser
from builder.databases.parsers.sider_parser import SIDERParser
from builder.databases.parsers.psp_parser import PhosphoSitePlusParser
from builder.databases.parsers.signor_parser import SIGNORParser
from builder.databases.parsers.corum_parser import CORUMParser
from builder.databases.parsers.exposome_parser import ExposomeExplorerParser
from builder.databases.parsers.hpa_parser import HPAParser


parsers = {
    # Base Databases
    "DrugBank": DrugBankParser,
    "UniProt": UniProtParser,
    "HGNC": HGNCParser,
    "HMDB": HMDBParser,
    "Reactome": ReactomeParser,
    "FooDB": FooDBParser,
    # DO is located on ontology part.
    # Domain Databases
    "CancerGenomeInterpreter": CancerGenomeInterpreterParser,
    "PathwayCommons": PathwayCommonsParser,
    "OncoKB": OncoKBParser,
    "Pfam": PfamParser,
    "RefSeq": RefSeqParser,
    "DisGEnet": DisGEnetParser,
    "STITCH": STITCHParser,
    "STRING": STRINGParser,
    "IntAct": IntActParser,
    "MutationDs": MutationDsParser,
    "JensenLab": JensenLabParser,
    "GWASCatalog": GWASCatalogParser,
    "DGIdb": DGIdbParser,
    "SIDER": SIDERParser,
    "PhosphoSitePlus": PhosphoSitePlusParser,
    "SIGNOR": SIGNORParser,
    "CORUM": CORUMParser,
    "ExposomeExplorer": ExposomeExplorerParser,
    "HPA": HPAParser
}
