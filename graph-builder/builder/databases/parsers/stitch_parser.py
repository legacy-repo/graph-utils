import os.path
import gzip
import csv
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class STITCHParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'STITCH'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self, drug_source=None):
        mapping = self.get_string_mapping()
        stored = set()
        relationship = None
        cutoff = self.config['STITCH_cutoff']
        header = self.config['header']
        drugmapping = {}
        evidences = ["experimental", "prediction",
                     "database", "textmining", "score"]
        relationship = "COMPILED_INTERACTS_WITH"
        url = self.config['STITCH_url']
        outputfile = os.path.join(
            self.import_directory, "stitch_associated_with.tsv")

        drugmapping = self.get_string_mapping(source=drug_source)
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])

        if self.download:
            self.download_db(url, directory)

        f = os.path.join(directory, fileName)
        associations = gzip.open(f, 'r')
        first = True
        with open(outputfile, 'w') as csvfile:
            writer = csv.writer(
                csvfile, delimiter='\t', escapechar='\\', quotechar='"', quoting=csv.QUOTE_ALL)
            writer.writerow(header)
            for line in associations:
                if first:
                    first = False
                    continue
                data = line.decode('utf-8').rstrip("\r\n").split()
                intA = data[0]
                intB = data[1]
                scores = data[2:]
                fscores = [str(float(score)/1000) for score in scores]

                if intA in drugmapping and intB in mapping and float(fscores[-1]) >= cutoff:
                    for aliasA in drugmapping[intA]:
                        for aliasB in mapping[intB]:
                            if (aliasA, aliasB) not in stored:
                                row = (aliasA, aliasB, relationship, "association",
                                       self.database_name, ",".join(evidences), ",".join(fscores[0:-1]), fscores[-1])
                                stored.add((aliasA, aliasB))
                                stored.add((aliasB, aliasB))
                                writer.writerow(row)
        associations.close()

        return mapping, drugmapping

    def parse_actions(self, proteinMapping, drugMapping=None):
        url = None
        header = self.config['header_actions']
        relationship = "COMPILED_ACTS_ON"
        stored = set()
        url = self.config['STITCH_actions_url']
        outputfile = os.path.join(self.import_directory,
                                  "stitch_drug_acts_on_protein.tsv")

        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)

        f = os.path.join(directory, fileName)
        associations = gzip.open(f, 'r')
        first = True
        with open(outputfile, 'w') as csvfile:
            writer = csv.writer(
                csvfile, delimiter='\t', escapechar='\\', quotechar='"', quoting=csv.QUOTE_ALL)
            writer.writerow(header)
            for line in associations:
                if first:
                    first = False
                    continue
                data = line.decode('utf-8').rstrip("\r\n").split()
                intA = data[0]
                intB = data[1]
                action = data[2]
                score = float(data[-1])/1000
                directionality = True

                if intB in proteinMapping:
                    aliasesA = []
                    if intA in drugMapping:
                        aliasesA = drugMapping[intA]
                    elif intA in proteinMapping:
                        aliasesA = proteinMapping[intA]
                    for aliasA in aliasesA:
                        for aliasB in proteinMapping[intB]:
                            if (aliasA, aliasB, action) not in stored:
                                row = (aliasA, aliasB, relationship,
                                       action, directionality, score, self.database_name)
                                writer.writerow(row)
                                stored.add((aliasA, aliasB, action))
                                stored.add((aliasB, aliasA, action))
        associations.close()

    def build_stats(self):
        proteinMapping, drugMapping = self.parse(
            drug_source=self.builder_config["database"]["sources"]["Drug"])
        self.parse_actions(proteinMapping, drugMapping)
        logger.success("Done Parsing database {}".format(self.database_name))
        return set()

    def get_string_mapping(self, source="BLAST_UniProt_AC"):
        """
        Parses database (db) and extracts relationships between identifiers to order databases (source).

        :param str url: link to download database raw file.
        :param str source: name of the source database for selecting aliases.
        :param bool download: wether to download the file or not.
        :param str db: name of the database to be parsed.
        :return: Dictionary of database identifers (keys) and set of unique aliases to other databases (values).
        """
        url = self.config['STITCH_mapping_url']
        mapping = defaultdict(set)
        directory = os.path.join(self.database_directory, self.database_name)
        file_name = os.path.join(directory, url.split('/')[-1])
        self.check_directory(directory)
        if self.download:
            logger.info("Downloading %s to %s" % (url, directory))
            self.download_db(url, directory)

        f = os.path.join(directory, file_name)
        first = True
        with gzip.open(f, 'rb') as mf:
            for line in mf:
                if first:
                    first = False
                    continue
                data = line.decode('utf-8').rstrip("\r\n").split("\t")
                stringID = data[0]
                alias = data[2]
                sources = data[3].split(' ')
                if not alias.startswith('DB'):
                    continue

                if source in sources:
                    mapping[stringID].add(alias)

        return mapping
