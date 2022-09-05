import os.path
import gzip
import re
import logging
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser
from builder.databases.parsers.stitch_parser import STITCHParser

logger = logging.getLogger(__name__)


class SIDERParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'SIDER'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self, drug_source):
        url = self.config['SIDER_url']
        header = self.config['header']

        output_file = 'sider_has_side_effect.tsv'

        stitch_parser = STITCHParser(self.import_directory, self.database_directory,
                                     download=self.download, skip=self.skip)
        drugmapping = stitch_parser.get_string_mapping(source=drug_source)
        phenotypemapping = self.get_mapping_from_ontology(ontology="Phenotype",
                                                          source=self.config['SIDER_source'])

        relationships = set()
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)
        associations = gzip.open(fileName, 'r')
        for line in associations:
            data = line.decode('utf-8').rstrip("\r\n").split("\t")
            drug = re.sub(r'CID\d', 'CIDm', data[0])
            se = data[2]
            evidence_from = str(data[3])
            #freq = data[4]
            #lower_bound = data[5]
            #upper_bound = data[6]
            if se.lower() in phenotypemapping and drug in drugmapping:
                for d in drugmapping[drug]:
                    p = phenotypemapping[se.lower()]
                    relationships.add(
                        (d, p, "HAS_SIDE_EFFECT", "SIDER", se, evidence_from))
        associations.close()

        return (relationships, header, output_file, drugmapping, phenotypemapping)

    def parse_indications(self, drugMapping, phenotypeMapping):
        url = self.config['SIDER_indications']
        header = self.config['indications_header']
        output_file = 'sider_is_indicated_for.tsv'

        relationships = set()
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)
        associations = gzip.open(fileName, 'r')
        for line in associations:
            data = line.decode('utf-8').rstrip("\r\n").split("\t")
            drug = re.sub(r'CID\d', 'CIDm', data[0])
            se = data[1]
            evidence = data[2]
            if se.lower() in phenotypeMapping and drug in drugMapping:
                for d in drugMapping[drug]:
                    p = phenotypeMapping[se.lower()]
                    relationships.add(
                        (d, p, "IS_INDICATED_FOR", evidence, "SIDER", se))

        associations.close()

        # self.remove_directory(directory)

        return (relationships, header, output_file)

    def build_stats(self):
        stats = set()
        source = self.builder_config["database"]["sources"]["Drug"]
        relationships,header, outputfileName, drugMapping, phenotypeMapping = self.parse(source)
        outputfile = os.path.join(self.import_directory, outputfileName)
        self.write_relationships(relationships, header, outputfile)
        logger.info("Database {} - Number of {} relationships: {}".format(self.database_name, 
                                                                          "has_side_effect", len(relationships)))
        stats.add(self._build_stats(len(relationships), "relationships", "has_side_effect", 
                                    self.database_name, outputfile, self.updated_on))
        relationships, header, outputfileName = self.parse_indications(drugMapping, phenotypeMapping)
        outputfile = os.path.join(self.import_directory, outputfileName)
        self.write_relationships(relationships, header, outputfile)
        logger.info("Database {} - Number of {} relationships: {}".format(self.database_name, 
                                                                          "indicated_for", len(relationships)))
        stats.add(self._build_stats(len(relationships), "relationships", "indicated_for", 
                                           self.database_name, outputfile, self.updated_on))
        logger.info("Done Parsing database {}".format(self.database_name))
        return stats
            