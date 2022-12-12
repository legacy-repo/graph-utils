# OncoKB database
import os.path
import re
import logging
import verboselogs
from collections import defaultdict
from builder.databases import config
from builder.databases.parsers.base_parser import BaseParser


logger = verboselogs.VerboseLogger('root')


class OncoKBParser(BaseParser):
    def __init__(self, import_directory, database_directory, config_file=None,
                 download=True, skip=True, organisms=["9606", "10090"]) -> None:
        self.database_name = 'OncoKB'
        config_dir = os.path.dirname(os.path.abspath(config.__file__))
        self.config_fpath = os.path.join(
            config_dir, "%s.yml" % self.database_name)

        super().__init__(import_directory, database_directory,
                         config_file, download, skip, organisms)

    def parse(self):
        url_actionable = self.config['OncoKB_actionable_url']
        url_annotation = self.config['OncoKB_annotated_url']
        amino_acids = self.config['amino_acids']
        entities_header = self.config['entities_header']
        relationships_headers = self.config['relationships_headers']
        mapping = self.get_mapping_from_ontology(
            ontology="Disease", source=None)

        drug_mapping = self.get_mapping_for_entity("Drug")
        protein_mapping = self.get_multiple_mapping_for_entity("Protein")

        levels = self.config['OncoKB_levels']
        entities = set()
        relationships = defaultdict(set)
        directory = os.path.join(self.database_directory, "OncoKB")
        self.check_directory(directory)
        acfileName = os.path.join(directory, os.path.basename(
            os.path.normpath(url_actionable)))
        anfileName = os.path.join(directory, os.path.basename(
            os.path.normpath(url_annotation)))

        if self.download:
            self.download_db(url_actionable, directory)
            self.download_db(url_annotation, directory)

        variant_regex = r"(\D\d+\D)$"
        with open(anfileName, 'r', errors='replace') as variants:
            first = True
            for line in variants:
                if line.startswith('#'):
                    continue
                elif first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                if len(data) >= 7:
                    gene = data[3]
                    variant = data[5]
                    oncogenicity = data[6]
                    effect = data[7]
                    if gene in protein_mapping:
                        for protein in protein_mapping[gene]:
                            match = re.search(variant_regex, variant)
                            if match:
                                if variant[0] in amino_acids and variant[-1] in amino_acids:
                                    valid_variant = protein + '_p.' + \
                                        amino_acids[variant[0]] + ''.join(
                                            variant[1:-1]) + amino_acids[variant[-1]]
                                    entities.add(
                                        (valid_variant, "Clinically_relevant_variant", "", "", "", "", "", effect, oncogenicity))

        with open(acfileName, 'r', errors='replace') as associations:
            first = True
            for line in associations:
                if line.startswith('#'):
                    continue
                elif first:
                    first = False
                    continue
                data = line.rstrip("\r\n").split("\t")
                if len(data) >= 9:
                    isoform = data[1]
                    gene = data[3]
                    variant = data[5]
                    disease = data[6]
                    level = data[7]
                    drugs = data[8].split(', ')
                    pubmed_ids = data[9].split(',')
                    if level in levels:
                        level = levels[level]

                    valid_variants = []
                    if gene in protein_mapping:
                        for protein in protein_mapping[gene]:
                            match = re.search(variant_regex, variant)
                            if match:
                                if variant[0] in amino_acids and variant[-1] in amino_acids:
                                    valid_variants.append(
                                        protein + '_p.' + amino_acids[variant[0]] + ''.join(variant[1:-1]) + amino_acids[variant[-1]])
                    for drug in drugs:
                        for d in drug.split(' + '):
                            if d.lower() in drug_mapping:
                                drug = drug_mapping[d.lower()]
                                relationships["targets"].add(
                                    (drug, gene, "CURATED_TARGETS", "curated", "NA", "NA", "curated", "OncoKB"))
                                for valid_variant in valid_variants:
                                    relationships["targets_clinically_relevant_variant"].add(
                                        (drug, valid_variant, "TARGETS_KNOWN_VARIANT", level[0], level[1], disease, "curated", "OncoKB"))
                    for valid_variant in valid_variants:
                        if disease.lower() in mapping:
                            disease = mapping[disease.lower()]
                            relationships["associated_with"].add(
                                (valid_variant, disease, "ASSOCIATED_WITH", "curated", "curated", "OncoKB", len(pubmed_ids)))
                        else:
                            pass
                        relationships["known_variant_is_clinically_relevant"].add(
                            (valid_variant, valid_variant, "KNOWN_VARIANT_IS_CLINICALLY_RELEVANT", "OncoKB"))

        # self.remove_directory(directory)

        return (entities, relationships, entities_header, relationships_headers)

    def build_stats(self):
        stats = set()
        entities, relationships, entities_header,  relationships_headers = self.parse()
        outputfile = os.path.join(
            self.import_directory, "oncokb_Clinically_relevant_variant.tsv")
        self.write_entities(entities, entities_header, outputfile)
        logger.info("Database {} - Number of {} entities: {}".format(
            self.database_name, "Clinically_relevant_variant", len(entities)))
        stats.add(self._build_stats(len(entities), "entity", "Clinically_relevant_variant",
                  self.database_name, outputfile, self.updated_on))
        for relationship in relationships:
            oncokb_outputfile = os.path.join(
                self.import_directory, "oncokb_"+relationship+".tsv")
            if relationship in relationships_headers:
                header = relationships_headers[relationship]
            else:
                header = ['START_ID', 'END_ID', 'TYPE']
            self.write_relationships(
                relationships[relationship], header, oncokb_outputfile)
            logger.info("Database {} - Number of {} relationships: {}".format(
                self.database_name, relationship, len(relationships[relationship])))
            stats.add(self._build_stats(len(relationships[relationship]), "relationships",
                                        relationship, self.database_name, outputfile, self.updated_on))
        logger.success("Done Parsing database {}".format(self.database_name))
        return stats
