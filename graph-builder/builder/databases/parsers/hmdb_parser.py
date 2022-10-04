# Human Metabolome Database
import os.path
import logging
import zipfile
import verboselogs
from collections import defaultdict
from lxml import etree
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class HMDBParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None, download=True, skip=True) -> None:
        self.database_name = 'HMDB'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory, config_file, download, skip)

    def parse(self):
        directory = os.path.join(self.database_directory, "HMDB")
        self.check_directory(directory)
        metabolites = self.extract_metabolites(directory)
        mapping = self.get_mapping_from_ontology(
            ontology="Disease", source=self.config['HMDB_DO_source'])
        mapping.update(self.get_mapping_from_ontology(
            ontology="Tissue", source=None))
        entities, attributes = self.build_metabolite_entity(directory, metabolites)
        relationships = self.build_relationships_from_hmdb(metabolites, mapping)
        entities_header = ['ID'] + attributes
        relationships_header = self.config['relationships_header']

        # self.remove_directory(directory)

        return (entities, relationships, entities_header, relationships_header)

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_header = self.parse()
        entity_outputfile = os.path.join(
            self.import_directory, "Metabolite.tsv")
        self.write_entities(entities, entities_header, entity_outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Metabolite", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Metabolite",
                  self.database_name, entity_outputfile, self.updated_on))
        for relationship in relationships:
            hmdb_outputfile = os.path.join(
                self.import_directory, relationship+".tsv")
            self.write_relationships(
                relationships[relationship], relationships_header, hmdb_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, relationship, len(relationships[relationship])))
            stats.add(self._build_stats(len(relationships[relationship]), "relationships",
                      relationship, self.database_name, hmdb_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats

    def extract_metabolites(self, directory):
        metabolites = defaultdict()
        prefix = "{http://www.hmdb.ca}"
        url = self.config['HMDB_url']
        fileName = os.path.join(directory, url.split('/')[-1])
        if self.download:
            self.download_db(url, directory)
        fields = self.config['HMDB_fields']
        parentFields = self.config['HMDB_parentFields']
        structuredFields = self.config['HMDB_structures']
        with zipfile.ZipFile(fileName, 'r') as zipped:
            for zfile in zipped.namelist():
                zipped.extract(member=zfile, path=directory)
                xfile = os.path.join(directory, zfile)
                with open(xfile, 'rb') as f:
                    context = etree.iterparse(f, events=(
                        "end",), tag=prefix + "metabolite")
                    for _, elem in context:
                        values = {child.tag.replace(prefix, ''): child.text for child in elem.iterchildren(
                        ) if child.tag.replace(prefix, '') in fields and child.text is not None}
                        for child in elem.iterchildren():
                            if child.tag.replace(prefix, '') in parentFields:
                                label = child.tag.replace(prefix, '')
                                values[label] = set()
                                for intchild in child.iter():
                                    if intchild.text is not None:
                                        text = intchild.text
                                        if text.strip() != "":
                                            if label in structuredFields:
                                                if intchild.tag.replace(prefix, '') in structuredFields[label]:
                                                    if len(structuredFields[label]) > 1:
                                                        values[intchild.tag.replace(
                                                            prefix, '')] = text
                                                    else:
                                                        values[label].add(text)
                                            elif intchild.tag.replace(prefix, '') in fields and text:
                                                values[label].add(text)
                        if "accession" in values:
                            metabolites[values["accession"]] = values

        return metabolites

    def build_metabolite_entity(self, directory, metabolites):
        entities = set()
        attributes = self.config['HMDB_attributes']
        for metid in metabolites:
            entity = []
            entity.append(metid)
            for attr in attributes:
                if attr in metabolites[metid]:
                    if type(metabolites[metid][attr]) == set:
                        lattr = ";".join(list(metabolites[metid][attr]))
                        entity.append(lattr)
                    else:
                        entity.append(metabolites[metid][attr])
                else:
                    entity.append('')
            entities.add(tuple(entity))

        self.build_hmdb_dictionary(directory, metabolites)

        return entities, attributes

    def build_relationships_from_hmdb(self, metabolites, mapping):
        relationships = defaultdict(list)
        associations = self.config['HMDB_associations']
        for metid in metabolites:
            for ass in associations:
                ident = ass
                if len(associations[ass]) > 1:
                    ident = associations[ass][1]
                if ass in metabolites[metid]:
                    if type(metabolites[metid][ass]) == set:
                        for partner in metabolites[metid][ass]:
                            if partner.lower() in mapping:
                                partner = mapping[partner.lower()]
                            relationships[ident].append(
                                (metid, partner, associations[ass][0], "HMDB"))
                    else:
                        partner = metabolites[metid][ass]
                        if metabolites[metid][ass].lower() in mapping:
                            partner = mapping[metabolites[metid][ass].lower()]
                        relationships[ident].append(
                            (metid, partner, associations[ass][0], "HMDB"))

        return relationships

    def build_hmdb_dictionary(self, directory, metabolites):
        filename = "mapping.tsv"
        outputfile = os.path.join(directory, filename)
        self.reset_mapping()
        with open(outputfile, 'w', encoding='utf-8') as out:
            for metid in metabolites:
                if "name" in metabolites[metid]:
                    name = metabolites[metid]["name"]
                    out.write(metid+"\t"+name.lower()+"\n")
                if "synonyms" in metabolites[metid]:
                    for synonym in metabolites[metid]["synonyms"]:
                        out.write(metid+"\t"+synonym.lower()+"\n")
                if "chebi_id" in metabolites[metid]:
                    chebi_id = metabolites[metid]["chebi_id"]
                    out.write(metid+"\t"+chebi_id+"\n")

        self.mark_complete_mapping()
