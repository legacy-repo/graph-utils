import os
import logging
import verboselogs
import zipfile
from os import path as os_path
from collections import defaultdict
from lxml import etree
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class NoSuchFile(Exception):
    pass


class DrugBankParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'DrugBank'
        config_dir = os_path.dirname(os_path.abspath(config.__file__))
        self.config_fpath = os_path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        directory = os_path.join(self.database_directory, "DrugBank")
        self.check_directory(directory)
        drugs = self.extract_drugs(directory)
        self.build_drug_bank_dictionary(directory, drugs)
        relationships = self.build_relationships_from_drug_bank(drugs)
        entities, attributes = self.build_drug_entity(drugs)
        entities_header = ['ID'] + attributes
        relationships_headers = self.config['relationships_headers']

        return (entities, relationships, entities_header, relationships_headers)

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header, relationships_headers = self.parse()
        entity_outputfile = os_path.join(self.import_directory, "Drug.tsv")
        self.write_entities(
            entities, entities_header, entity_outputfile)
        logger.info(
            "Database {} - Number of {} entities: {}".format(self.database_name, "Drug", len(entities)))
        stats.add(self._build_stats(len(entities), "entity",
                  "Drug", self.database_name, entity_outputfile, self.updated_on))
        for relationship in relationships:
            relationship_outputfile = os_path.join(
                self.import_directory, relationship+".tsv")
            header = ['START_ID', 'END_ID', 'TYPE', 'source']
            if relationship in relationships_headers:
                header = relationships_headers[relationship]
            self.write_relationships(relationships[relationship],
                                     header, relationship_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, relationship, len(relationships[relationship])))
            stats.add(self._build_stats(len(relationships[relationship]),
                                        "relationships", relationship, self.database_name,
                                        relationship_outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats

    def extract_drugs(self, directory):
        drugs = {}
        prefix = '{http://www.drugbank.ca}'
        url = self.config['DrugBank_url']
        fileName = os_path.join(directory, url.split('/')[-1])
        fields = self.config['DrugBank_fields']
        attributes = self.config['DrugBank_attributes']
        parentFields = self.config['DrugBank_parentFields']
        structuredFields = self.config['DrugBank_structures']

        self.check_file(fileName)

        vocabulary = self.parse_drug_bank_vocabulary(directory)
        with zipfile.ZipFile(fileName, 'r') as zipped:
            for zfile in zipped.namelist():
                zipped.extract(member=zfile, path=directory)
                xfile = os_path.join(directory, zfile)
                with open(xfile, 'rb') as f:
                    context = etree.iterparse(
                        f, events=("end",), tag=prefix+"drug")
                    for a, elem in context:
                        synonyms = set()
                        values = {child.tag.replace(prefix, ''): child.text for child in elem.iterchildren(
                        ) if child.tag.replace(prefix, '') in fields and child.text is not None}
                        if "drugbank-id" in values:
                            synonyms.add(values["drugbank-id"])
                        for child in elem.iterchildren():
                            if child.tag.replace(prefix, '') in parentFields:
                                label = child.tag.replace(prefix, '')
                                values[label] = []
                                for intchild in child.iter():
                                    if intchild.text is not None and intchild.text.strip() != "":
                                        if label in structuredFields:
                                            if intchild.tag.replace(prefix, '') in structuredFields[label]:
                                                if label == "external-identifiers":
                                                    synonyms.add(intchild.text)
                                                else:
                                                    values[label].append(
                                                        intchild.text)
                                        elif intchild.tag.replace(prefix, '') in fields and intchild.text:
                                            values[label].append(intchild.text)
                                        elif intchild.tag.replace(prefix, '') in attributes and intchild.text:
                                            values[intchild.tag.replace(
                                                prefix, '')] = intchild.text

                        if "drugbank-id" in values and len(values) > 2:
                            if values["drugbank-id"] in vocabulary:
                                values["id"] = vocabulary[values["drugbank-id"]]
                                synonyms.add(values["drugbank-id"])
                                # values["alt_drugbank-id"] = vocabulary[values['id']]
                                values["synonyms"] = list(synonyms)
                                drugs[values["id"]] = values
        return drugs

    def check_file(self, filepath):
        if not os.path.exists(filepath):
            raise NoSuchFile(
                "Cannot locate the %s, please download it manually(need to login)." % filepath)

    def parse_drug_bank_vocabulary(self, directory):
        vocabulary = {}
        url = self.config['DrugBank_vocabulary_url']
        fileName = os_path.join(directory, url.split('/')[-1])

        self.check_file(fileName)
        with zipfile.ZipFile(fileName, 'r') as zipped:
            for f in zipped.namelist():
                with zipped.open(f) as vf:
                    # with open(os.path.join(directory,f), 'r') as vf:
                    for line in vf:
                        data = line.decode('utf-8').rstrip('\r\n').split(',')
                        primary = data[0]
                        secondaries = data[1].split(' | ')
                        for sec in secondaries:
                            vocabulary[sec] = primary
                            vocabulary[primary] = primary
        return vocabulary

    def build_relationships_from_drug_bank(self, drugs):
        relationships = defaultdict(list)
        associations = self.config['DrugBank_associations']
        for did in drugs:
            for ass in associations:
                ident = ass
                if len(associations[ass]) > 1:
                    ident = associations[ass][1]
                if ass in drugs[did]:
                    if type(drugs[did][ass]) == list:
                        partners = drugs[did][ass]
                        if ass == "drug-interactions":
                            partners = zip(partners[0::2], partners[1::2])
                        elif ass in ["snp-effects", 'snp-adverse-drug-reactions']:
                            partners = zip(
                                partners[0::3], partners[1::3], partners[2::3])
                        elif ass == "targets":
                            partners = zip(partners[0::2], partners[1::2])
                            partners = [
                                p for r, p in partners if r == "UniProtKB"]
                        for partner in partners:
                            rel = (did, partner,
                                   associations[ass][0], "DrugBank")
                            relationships[ident].append(
                                tuple(self.flatten(rel)))
                    else:
                        partner = drugs[did][ass]
                        relationships[ident].append(
                            (did, partner, associations[ass][0], "DrugBank"))

        return relationships

    def build_drug_entity(self, drugs):
        entities = set()
        attributes = self.config['DrugBank_attributes']
        properties = self.config['DrugBank_exp_prop']
        allAttr = attributes[:-1] + [p.replace(' ', '_').replace(
            '-', '_').replace('(', '').replace(')', '') for p in properties]
        for did in drugs:
            entity = []
            entity.append(did)
            for attr in attributes:
                if attr in drugs[did]:
                    if type(drugs[did][attr]) == list:
                        if attr == "experimental-properties":
                            newAttr = dict(
                                zip(drugs[did][attr][0::2], drugs[did][attr][1::2]))
                            for prop in properties:
                                if prop in newAttr:
                                    entity.append(newAttr[prop])
                                else:
                                    entity.append('')
                        else:
                            lattr = "|".join(drugs[did][attr])
                            entity.append(lattr)
                    else:
                        entity.append(drugs[did][attr])
                else:
                    entity.append('')
            entities.add(tuple(entity))

        return entities, allAttr

    def build_drug_bank_dictionary(self, directory, drugs):
        filename = self.config['DrugBank_dictionary_file']
        outputfile = os_path.join(directory, filename)
        self.reset_mapping()
        with open(outputfile, 'w', encoding='utf-8') as out:
            for did in drugs:
                if "name" in drugs[did]:
                    name = drugs[did]["name"]
                    out.write(did+"\t"+name.lower()+"\n")
                if "synonyms" in drugs[did]:
                    for synonym in drugs[did]["synonyms"]:
                        out.write(did+"\t"+synonym.lower()+"\n")

        self.mark_complete_mapping()
