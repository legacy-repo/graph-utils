import os.path
import logging
import verboselogs
import gzip
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class DisGEnetParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'DisGEnet'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        relationships = defaultdict(set)

        files = self.config['disgenet_files']
        mapping_files = self.config['disgenet_mapping_files']
        url = self.config['disgenet_url']
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        header = self.config['disgenet_header']
        output_file = 'disgenet_associated_with.tsv'

        if self.download:
            for f in files:
                self.download_db(url+files[f], directory)
            for f in mapping_files:
                self.download_db(url+mapping_files[f], directory)

        protein_mapping = self.read_disgenet_protein_mapping(directory)
        disease_mapping = self.read_disgenet_disease_mapping(directory)
        for f in files:
            first = True
            associations = gzip.open(os.path.join(directory, files[f]), 'r')
            dtype, atype = f.split('_')
            if dtype == 'gene':
                idType = "Protein"
                scorePos = 9
            if dtype == 'variant':
                idType = "Transcript"
                scorePos = 5
            for line in associations:
                if first:
                    first = False
                    continue
                try:
                    data = line.decode('utf-8').rstrip("\r\n").split("\t")
                    geneId = str(int(data[0]))
                    #disease_specificity_index =  data[2]
                    #disease_pleiotropy_index = data[3]
                    diseaseId = data[4]
                    score = float(data[scorePos])
                    pmids = data[13]
                    source = data[-1]
                    if geneId in protein_mapping:
                        for identifier in protein_mapping[geneId]:
                            if diseaseId in disease_mapping:
                                for code in disease_mapping[diseaseId]:
                                    code = "DOID:"+code
                                    relationships[idType].add(
                                        (identifier, code, "ASSOCIATED_WITH", score, atype, "DisGeNet: "+source, pmids))
                except UnicodeDecodeError:
                    continue
            associations.close()

        # self.remove_directory(directory)

        return (relationships, header, output_file)

    def read_disgenet_protein_mapping(self, directory):
        files = self.config['disgenet_mapping_files']
        first = True
        mapping = defaultdict(set)
        if "protein_mapping" in files:
            mappingFile = files["protein_mapping"]
            with gzip.open(os.path.join(directory, mappingFile), 'r') as f:
                for line in f:
                    if first:
                        first = False
                        continue
                    data = line.decode('utf-8').rstrip("\r\n").split("\t")
                    identifier = data[0]
                    intIdentifier = data[1]
                    mapping[intIdentifier].add(identifier)
        return mapping

    def read_disgenet_disease_mapping(self, directory):
        files = self.config['disgenet_mapping_files']
        first = True
        mapping = defaultdict(set)
        if "disease_mapping" in files:
            mappingFile = files["disease_mapping"]
            with gzip.open(os.path.join(directory, mappingFile), 'r') as f:
                for line in f:
                    if first:
                        first = False
                        continue
                    data = line.decode('utf-8').rstrip("\r\n").split("\t")
                    identifier = data[0]
                    vocabulary = data[2]
                    code = data[3]
                    if vocabulary == "DO":
                        mapping[identifier].add(code)
        return mapping

    def build_stats(self):
        stats = set()
        relationships, header, outputfileName = self.parse()
        for idType in relationships:
            outputfile = os.path.join(self.import_directory,
                                      idType+"_"+outputfileName)
            self.write_relationships(relationships[idType], header, outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, idType, len(relationships[idType])))
            stats.add(self._build_stats(len(relationships[idType]), "relationships", idType,
                                        self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
