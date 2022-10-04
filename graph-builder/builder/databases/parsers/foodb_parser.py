# FooDB
import os.path
import tarfile
import pandas as pd
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class FooDBParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'FooDB'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, "FooDB")
        self.check_directory(directory)

        database_url = self.config['database_url']
        entities_header = self.config['entities_header']
        relationships_headers = self.config['relationships_headers']
        tar_fileName = os.path.join(directory, database_url.split('/')[-1])
        if self.download:
            self.download_db(database_url, directory)

        contents = {}
        food = set()
        compounds = {}
        try:
            tf = tarfile.open(tar_fileName, 'r')
            file_content = tf.getnames()
            tar_dir = file_content[1]
            tf.extractall(path=directory)
            tf.close()
            for file_name in self.config['files']:
                path = os.path.join(
                    directory, os.path.join(tar_dir, file_name))
                with open(path, 'r', encoding="utf-8", errors='replace') as f:
                    if file_name == "Content.csv":
                        contents = self.parse_contents(f)
                    elif file_name == "Food.csv":
                        food, mapping = self.parse_food(f)
                    elif file_name == "Compound.csv":
                        compounds = self.parse_compounds(f)
            for food_id, compound_id in contents:
                if compound_id in compounds:
                    compound_code = compounds[compound_id].replace(
                        "HMDB", "HMDB00")
                    relationships[("food", "has_content")].add(
                        (food_id, compound_code, "HAS_CONTENT") + contents[(food_id, compound_id)])
            self.reset_mapping()
            with open(os.path.join(directory, "mapping.tsv"), 'w', encoding='utf-8') as out:
                for food_id in mapping:
                    for alias in mapping[food_id]:
                        out.write(str(food_id)+"\t"+str(alias)+"\n")

            self.mark_complete_mapping()
        except tarfile.ReadError as err:
            raise Exception("Error importing database FooDB.\n {}".format(err))

        # self.remove_directory(directory)

        return food, relationships, entities_header, relationships_headers

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_headers = self.parse()
        entity_outputfile = os.path.join(self.import_directory, "Food.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Food", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Food",
                  self.database_name, entity_outputfile, self.updated_on))
        for entity, relationship in relationships:
            foodb_outputfile = os.path.join(self.import_directory,
                                            self.database_name.lower()+"_"+entity.lower()+"_"+relationship.lower()+".tsv")
            self.write_relationships(relationships[(entity, relationship)],
                                     relationships_headers[entity], foodb_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, relationship, len(relationships[(entity, relationship)])))
            stats.add(self._build_stats(len(relationships[(entity, relationship)]),
                                        "relationships", relationship, self.database_name, foodb_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats

    def parse_contents(self, fhandler):
        contents = {}
        first = True
        for line in fhandler:
            if first:
                first = False
                continue
            data = line.rstrip("\r\n").split(",")
            if len(data) == 24:
                compound_id = data[0]
                food_id = int(data[3])
                min_cont = float(data[11]) if data[11] != 'NULL' else 0
                max_cont = float(data[12]) if data[12] != 'NULL' else 0
                units = data[13].replace('"', '')
                average = float(data[23]) if data[23] != 'NULL' else 0
                contents[(food_id, compound_id)] = (
                    min_cont, max_cont, average, units, "FooDB")
        return contents

    def parse_food(self, fhandler):
        food = set()
        mapping = defaultdict(set)
        df = pd.read_csv(fhandler, sep=',', header=None,
                         low_memory=False, encoding="utf-8")
        first = True
        for index, row in df.iterrows():
            if first:
                first = False
                continue
            food_id = row[22]
            name = row[1]
            sci_name = row[2]
            description = str(row[3]).replace('"', '')
            group = row[11]
            subgroup = row[12]
            food.add((food_id, name, sci_name, description,
                     group, subgroup, "FooDB"))
            mapping[food_id].add(name)
            mapping[food_id].add(sci_name)

        return food, mapping

    def parse_compounds(self, fhandler):
        compounds = {}
        first = True
        df = pd.read_csv(fhandler, sep=',', header=None,
                         low_memory=False, encoding="utf-8")
        first = True
        for index, row in df.iterrows():
            if first:
                first = False
                continue
            print(row)
            print(row.shape)
            compound_id = row[0]
            mapped_code = row[44]
            if str(mapped_code) != 'nan':
                compounds[compound_id] = mapped_code
        return compounds
