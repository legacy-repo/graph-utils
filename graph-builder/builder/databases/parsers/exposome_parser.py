import os.path
import zipfile
import pandas as pd
import os.path
import zipfile
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class ExposomeExplorerParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'ExposomeExplorer'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        directory = os.path.join(self.database_directory, self.database_name)
        self.check_directory(directory)
        database_urls = self.config['database_urls']
        relationships_header = self.config['relationships_header']
        mapping = self.get_mapping_for_entity("Food")
        correlations = {}
        for url in database_urls:
            zipped_fileName = os.path.join(directory, url.split('/')[-1])
            file_name = '.'.join(url.split('/')[-1].split('.')[0:2])
            if self.download:
                self.download_db(url, directory)

            with zipfile.ZipFile(zipped_fileName) as z:
                if file_name == "biomarkers.csv":
                    biomarkers = self.parse_biomarkers_file(z, file_name)
                elif file_name == "correlations.csv":
                    correlations = self.parse_correlations_file(z, file_name,
                                                                biomarkers, mapping)
        # self.remove_directory(directory)

        return correlations, relationships_header

    def parse_biomarkers_file(self, fhandler, file_name):
        biomarkers = {}
        first = True
        with fhandler.open(file_name) as f:
            df = pd.read_csv(f, sep=',', header=None,
                             low_memory=False)
            first = True
            for index, row in df.iterrows():
                if first:
                    first = False
                    continue
                identifier = row[0]
                metabolite = row[10]
                if metabolite != '':
                    biomarkers[identifier] = metabolite

        return biomarkers

    def parse_correlations_file(self, fhandler, file_name, biomarkers, mapping):
        correlations = defaultdict(set)
        first = True
        with fhandler.open(file_name) as f:
            df = pd.read_csv(f, sep=',', header=None,
                             low_memory=False)
            first = True
            for index, row in df.iterrows():
                if first:
                    first = False
                    continue
                biomarker = row[0]
                food_name = row[9]
                intake_median = row[14]
                intake_units = row[15]
                biosample = row[18]
                method = row[19]
                #corr_method = row[29]
                corr = float(row[29])
                ci_low = row[30]
                ci_high = row[31]
                pvalue = row[32]
                significant = row[33]
                publication = row[37]

                if significant in ["Yes", "yes", "YES", "Y"]:
                    if food_name in mapping:
                        food_id = mapping[food_name]
                        if biomarker in biomarkers:
                            biomarker_id = biomarkers[biomarker]
                            correlations[("food", "correlated_with_metabolite")].add((food_id, biomarker_id, "CORRELATED_WITH_METABOLITE", intake_median,
                                                                                      intake_units, biosample, method, corr, ci_low, ci_high, pvalue, significant, publication, "Exposome Explorer"))

        return correlations

    def build_stats(self):
        stats = set()
        relationships, header = self.parse()
        for entity, relationship in relationships:
            ee_outputfile = os.path.join(self.import_directory,
                                         self.database_name.lower()+"_"+entity.lower()+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     header[entity], ee_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(self.database_name, relationship,
                                                                              len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]), "relationships",
                                        relationship, self.database_name, ee_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
