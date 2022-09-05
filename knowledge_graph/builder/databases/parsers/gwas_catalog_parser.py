import os.path
import re
from collections import defaultdict
import os.path
import logging
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser

logger = logging.getLogger(__name__)


class GWASCatalogParser(BaseParser):
    def __init__(self, import_directory, database_directory, download=True, skip=True) -> None:
        self.database_name = 'GWASCatalog'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, download, skip)

    def parse(self):
        url = self.config['GWASCat_url']
        entities_header = self.config['entities_header']
        relationships_header = self.config['relationships_header']
        entities = set()
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)
        with open(fileName, 'r', encoding="utf-8") as catalog:
            for line in catalog:
                data = line.rstrip("\r\n").split("\t")
                if len(data) > 36:
                    pubmedid = data[1]
                    date = data[3]
                    title = data[6]
                    sample_size = data[8]
                    replication_size = data[9]
                    #chromosome = data[11]
                    #position = data[12]
                    #genes_mapped = data[14].split(" - ")
                    snp_id = data[20].split('-')[0]
                    freq = data[26]
                    pval = data[27]
                    odds_ratio = data[30]
                    trait = data[34]
                    exp_factor = data[35]
                    study = data[36]

                    entities.add((study, "GWAS_study", title, date,
                                 sample_size, replication_size, trait))
                    if pubmedid != "":
                        relationships["published_in_publication"].add((study, pubmedid,
                                                                       "PUBLISHED_IN", "GWAS Catalog"))
                    if snp_id != "":
                        relationships["variant_found_in_gwas"].add((re.sub(r"^\W+|\W+$", "", snp_id), study,
                                                                    "VARIANT_FOUND_IN_GWAS", freq, pval, odds_ratio, trait, "GWAS Catalog"))
                    if exp_factor != "":
                        exp_factor = exp_factor.split('/')[-1]
                        exp_factor = exp_factor.replace('_', ':')
                        relationships["studies_trait"].add((study, exp_factor,
                                                            "STUDIES_TRAIT", "GWAS Catalog"))

        # self.remove_directory(directory)

        return (entities, relationships, entities_header, relationships_header)

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_header = self.parse()
        entity_outputfile = os.path.join(self.import_directory,
                                         "GWAS_study.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(self.database_name,
                                                                     "GWAS_study", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "GWAS_study", self.database_name,
                                    entity_outputfile, self.updated_on))
        for relationship in relationships:
            header = ['START_ID', 'END_ID', 'TYPE', 'source']
            if relationship in relationships_header:
                header = relationships_header[relationship]
            outputfile = os.path.join(self.import_directory,
                                      "GWAS_study_"+relationship+".tsv")
            self.write_relationships(relationships[relationship],
                                     header, outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(self.database_name,
                                                                              relationship, len(relationships[relationship])))
            stats.add(self._build_stats(len(relationships[relationship]), "relationships", relationship,
                                        self.database_name, outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
